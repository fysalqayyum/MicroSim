#ifndef FARFIELD_GUARD_H
#define FARFIELD_GUARD_H

/*
 * Runtime validity guards for directional-solidification runs with a moving
 * window.  Optional: nothing here changes the solution, it only observes it and
 * aborts on states that are known to produce physically meaningless output.
 *
 * WHY THIS EXISTS.  A directional-solidification run can finish "cleanly" by
 * every conventional check and still be garbage.  The motivating failure ran to
 * completion with max|sum(phi)-1| ~ 1e-13, zero non-finite values, and a moving
 * window mass audit exact to ~1e-13 over several hundred shift events -- while
 * the front had banked undercooling well past the freezing range and then
 * released it at more than ten times the pulling velocity into a grid-scale
 * structure.
 *
 * The lesson generalises: sum(phi) = 1 is satisfied by garbage, and a runaway is
 * perfectly finite.  A metric that cannot represent the failure mode is not a
 * guard.  The checks below watch quantities that CAN represent it.
 *
 * LIQUIDUS MODEL.  A linearised binary liquidus,
 *   T_liq(c) = T_liq(c0) + m_L * (wt%(c) - wt%(c0))
 * used only for reporting.  Every constant is read from the environment rather
 * than the .in file, so alloy data can be changed without touching the parser --
 * silently unparsed .in keys are a known failure mode.
 *
 *   FARFIELD_T_LIQ_C0   liquidus at c0, K
 *   FARFIELD_M_L        liquidus slope, K per wt% solute
 *   FARFIELD_C0_MOL     reference composition, mole fraction
 *   FARFIELD_BAND       physical rows averaged at the top (default 8)
 *
 * WARNING: the compiled-in fallbacks for the first three are Al-Si values, and
 * farfield_wt_si() assumes a binary Al-Si liquid.  Set all three explicitly for
 * any other alloy, and replace farfield_wt_si() if the solvent is not Al.  The
 * fallbacks exist so a missing export degrades to a wrong LOG LINE rather than a
 * crash; they are not a default worth relying on.
 *
 * Exit code 17 has ONE trigger: solid cells in the top band, i.e. the front has
 * reached the top of the window and the moving window can no longer track.
 * Unconditional.
 *
 * A far-field undercooling abort was tried and REMOVED.  In a linear thermal
 * gradient the far-field undercooling is an algebraic identity in the tip
 * undercooling and the box height, so it carries no information the tip state
 * does not already contain, and it merely re-expresses a fixed-headroom sizing
 * rule in kelvin -- unpassable by construction in a deliberately short box.  The
 * quantity is still computed and logged; only the abort is gone.  See the note
 * at the removal site.
 *
 * ------------------------------------------------------------------------
 * SHORT-BOX GUARDS
 *
 * Undercooled bulk liquid is INERT in this model: h(phi) = 6phi^5 - 15phi^4 +
 * 10phi^3 for npha = 2, so h'(0) = h''(0) = 0 and the driving-force term
 * vanishes to second order at phi = 0; and addNoise returns immediately when
 * envelope = phi0*phi1 <= 0, so bulk liquid receives exactly zero stochastic
 * forcing.  Uniform liquid at phi = 0 is therefore an exact stationary state at
 * ANY undercooling and cannot nucleate.  This has been confirmed in practice:
 * liquid held more than 14 K below its liquidus for over a million steps showed
 * no nucleation.
 *
 * Consequently a box may be sized for SOLUTE CONTAINMENT rather than to reach
 * the liquidus isotherm, which is much cheaper.  These are then the checks that
 * matter:
 *
 *   SHORTBOX_DC_MAX_REL  abort when the top band's mean composition rises this
 *                        far ABOVE c0 (default 0.02 = 2%).  MANDATORY for a
 *                        short box: with a no-flux top a truncated boundary
 *                        layer traps rejected solute and the alloy silently
 *                        enriches, and the moving-window mass audit CANNOT see
 *                        it -- its `expected_after` is defined to include
 *                        refill at c0.  This is the only check that looks at
 *                        the quantity that actually drifts.
 *   SHORTBOX_TIP_MARGIN  abort when the tip comes within this many cells of the
 *                        top boundary (default 200).  Fires earlier and more
 *                        informatively than "solid reached the top band".
 *   SHORTBOX_VTIP_FACTOR abort when v_tip exceeds this multiple of the pulling
 *                        velocity (default 2.0) -- the runaway signature.
 *                        Measured on the ABSOLUTE front position: tip row plus
 *                        shift_OFFSET plus shift_position, matching what
 *                        shift.dat records.  Using shift_position alone is
 *                        WRONG: it is set once when a restart file is read and
 *                        never updated, so the check is blinded by the very
 *                        moving window it exists to see through.
 *   SHORTBOX_GUARD_OFF   set to 1 to log only, never abort (all three)
 *
 * Exit codes 21/22/23.  Each check logs every call whether or not it aborts.
 */

#define FARFIELD_EXIT_CODE 17
#define SHORTBOX_EXIT_DC   21
#define SHORTBOX_EXIT_TIP  22
#define SHORTBOX_EXIT_VTIP 23

static double farfield_env_double(const char *name, double fallback) {
  const char *raw = getenv(name);
  char *end;
  double value;

  if (raw == NULL || raw[0] == '\0') {
    return fallback;
  }
  value = strtod(raw, &end);
  if (end == raw) {
    if (rank == MASTER) {
      printf("FARFIELD GUARD: cannot parse %s='%s', using %g\n",
             name, raw, fallback);
    }
    return fallback;
  }
  return value;
}

/* wt% Si of a binary Al-Si liquid at mole fraction x. */
static double farfield_wt_si(double x) {
  const double M_SI = 28.0855;
  const double M_AL = 26.9815;
  return 100.0 * x * M_SI / (x * M_SI + (1.0 - x) * M_AL);
}

/*
 * Highest physical grid row containing solid, or -1 if the domain is all
 * liquid.  Scans downward from the top and stops at the first row that has
 * any, so on a normal run it touches only the empty liquid rows.
 */
static long farfield_tip_row(void) {
  long x, y, z, idx;
  long local_tip = -1, global_tip;

  for (y = mpiparam.rows_y - 2; y >= 1 && local_tip < 0; y--) {
    for (x = 1; x < mpiparam.rows_x - 1 && local_tip < 0; x++) {
      for (z = 1; z < mpiparam.rows_z - 1; z++) {
        idx = y + mpiparam.rows_y * (z + mpiparam.rows_z * x);
        if (gridinfomN[idx].phia[NUMPHASES - 1] < 0.5) {
          local_tip = y;
          break;
        }
      }
    }
  }
  MPI_Allreduce(&local_tip, &global_tip, 1, MPI_LONG, MPI_MAX, MPI_COMM_WORLD);
  return global_tip;
}

/*
 * Averages the top FARFIELD_BAND physical rows and compares their temperature
 * against the liquidus of their own composition.  Logs one line per call.
 * Aborts only if solid has reached the band at all, which means the front has
 * run out of domain and the moving window can no longer track it.
 *
 * Then runs the three short-box checks (composition containment, tip margin,
 * front velocity) documented at the top of this file.
 */
void farfield_guard_check(long absolute_step) {
  long x, y, z, idx;
  long y_hi = mpiparam.rows_y - 2;
  long y_lo;
  long band;
  long local_n = 0, global_n;
  double local_T = 0.0, local_c = 0.0, local_solid = 0.0;
  double global_T, global_c, global_solid;
  double T_liq_c0, m_L, c0_mol, c0_wt;
  double mean_T, mean_c, T_liq, undercooling;
  /* short-box guards */
  static long prev_step = -1;
  static long prev_front = 0;
  long tip_row, front_abs, tip_margin_cells;
  double dc_rel, dc_max_rel, vtip = 0.0, vtip_factor, V_pull;
  int shortbox_off, have_vtip = 0;

  /* FARFIELD_DT_MAX_K and FARFIELD_GUARD_OFF were removed together with the
   * undercooling abort (see the note further down).  Exporting either is now a
   * harmless no-op, so old launch scripts keep working unchanged. */
  T_liq_c0  = farfield_env_double("FARFIELD_T_LIQ_C0", 899.604509);
  m_L       = farfield_env_double("FARFIELD_M_L", -6.507);
  c0_mol    = farfield_env_double("FARFIELD_C0_MOL", 0.052953);
  band      = (long)farfield_env_double("FARFIELD_BAND", 8.0);

  if (band < 1) {
    band = 1;
  }
  y_lo = y_hi - (band - 1);
  if (y_lo < 1) {
    y_lo = 1;
  }

  for (x = 1; x < mpiparam.rows_x - 1; x++) {
    for (z = 1; z < mpiparam.rows_z - 1; z++) {
      for (y = y_lo; y <= y_hi; y++) {
        idx = y + mpiparam.rows_y * (z + mpiparam.rows_z * x);
        local_T += gridinfomN[idx].temperature;
        local_c += gridinfomN[idx].composition[0];
        /* phia[NUMPHASES-1] is the liquid fraction in this build. */
        if (gridinfomN[idx].phia[NUMPHASES - 1] < 0.5) {
          local_solid += 1.0;
        }
        local_n++;
      }
    }
  }

  MPI_Allreduce(&local_T, &global_T, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&local_c, &global_c, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&local_solid, &global_solid, 1, MPI_DOUBLE, MPI_SUM,
                MPI_COMM_WORLD);
  MPI_Allreduce(&local_n, &global_n, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);

  if (global_n == 0) {
    return;
  }

  mean_T = global_T / (double)global_n;
  mean_c = global_c / (double)global_n;
  c0_wt = farfield_wt_si(c0_mol);
  T_liq = T_liq_c0 + m_L * (farfield_wt_si(mean_c) - c0_wt);
  undercooling = T_liq - mean_T;

  /*
   * Short-box quantities are computed HERE, before the log write, so that
   * every diagnostic lands in the same CSV row even if a far-field check
   * aborts below.  The corresponding abort tests run at the end.
   */
  dc_max_rel       = farfield_env_double("SHORTBOX_DC_MAX_REL", 0.02);
  tip_margin_cells = (long)farfield_env_double("SHORTBOX_TIP_MARGIN", 200.0);
  vtip_factor      = farfield_env_double("SHORTBOX_VTIP_FACTOR", 2.0);
  shortbox_off     = (getenv("SHORTBOX_GUARD_OFF") != NULL
                      && getenv("SHORTBOX_GUARD_OFF")[0] == '1');

  V_pull  = temperature_gradientY.velocity;
  tip_row = farfield_tip_row();
  /* Absolute front position: the moving window must not be able to mask a
   * runaway, so measure in lab cells, not grid cells.
   *
   * NOTE.  shift_position is set once, when a restart file is read
   * (microsim_kks_opencl.c), and never changes during a run; the live
   * accumulator is shift_OFFSET (CL_Shift.h: shift_OFFSET += shift_cells).
   * Using shift_position alone freezes this value as soon as the moving window
   * pins tip_row at Shiftj: v_tip then reads exactly 0 for the rest of the run
   * and the SHORTBOX_VTIP_FACTOR check below can never fire -- precisely the
   * regime it exists to police.  shift.dat records
   * `shift_OFFSET + shift_position`, so match it exactly. */
  front_abs = tip_row + shift_OFFSET + shift_position;
  if (prev_step >= 0 && absolute_step > prev_step) {
    vtip = (double)(front_abs - prev_front) * deltay
           / ((double)(absolute_step - prev_step) * deltat);
    have_vtip = 1;
  }
  prev_step  = absolute_step;
  prev_front = front_abs;

  dc_rel = (c0_mol != 0.0) ? (mean_c - c0_mol) / c0_mol : 0.0;

  if (rank == MASTER) {
    FILE *log = fopen("DATA/farfield_guard.csv", "a+");
    if (log != NULL) {
      fseek(log, 0, SEEK_END);
      if (ftell(log) == 0) {
        fprintf(log, "step,rows_averaged,mean_T_K,mean_xSi,mean_wt_Si,"
                     "T_liquidus_K,far_field_undercooling_K,solid_cells_in_band,"
                     "tip_row,front_abs_cells,cells_below_top,"
                     "topband_dc_rel_c0,v_tip_m_s,v_tip_over_V\n");
      }
      fprintf(log, "%ld,%ld,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,"
                   "%ld,%ld,%ld,%.17g,%.17g,%.17g\n",
              absolute_step, global_n, mean_T, mean_c,
              farfield_wt_si(mean_c), T_liq, undercooling, global_solid,
              tip_row, front_abs, (mpiparam.rows_y - 2) - tip_row,
              dc_rel,
              have_vtip ? vtip : NAN,
              (have_vtip && V_pull != 0.0) ? vtip / V_pull : NAN);
      fclose(log);
    }
    printf("  far-field: T = %.3f K, xSi = %.6f (%.3f wt%%), "
           "T_liq = %.3f K, undercooling = %+.3f K (diagnostic only)\n",
           mean_T, mean_c, farfield_wt_si(mean_c), T_liq, undercooling);
    fflush(stdout);
  }

  if (global_solid > 0.0) {
    if (rank == MASTER) {
      printf("\nFARFIELD GUARD FAILED at step %ld: %g of %ld cells in the top "
             "%ld rows are SOLID.\n"
             "  The front has reached the top of the window; there is no liquid "
             "left to grow into and the moving window can no longer track.\n"
             "  Increase ny (thermal headroom) or lower Shiftj.\n",
             absolute_step, global_solid, global_n, band);
      fflush(stdout);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    exit(FARFIELD_EXIT_CODE);
  }

  /*
   * REMOVED: an `undercooling > FARFIELD_DT_MAX_K` abort that used to sit here.
   *
   * It tested no independent quantity.  In a linear thermal gradient the
   * far-field undercooling is an algebraic identity:
   *
   *   undercooling == dT_tip - cells_above_tip*G*dx + (T_liq_band - T_liq_c0)
   *
   * verified over 45 steady-state samples at a residual of 0.0000 K, sd
   * 0.0000 K -- exact to machine precision.  So `undercooling > dT_max` reduces
   * to `cells_above_tip < (dT_tip - dT_max)/(G*dx)`, i.e. a fixed-headroom box
   * sizing rule restated in kelvin.  For a box deliberately sized for solute
   * containment rather than to reach the liquidus isotherm it can never pass,
   * and it will abort otherwise valid runs partway through their transient.
   *
   * Nor could it detect what its message claimed.  Banking of undercooling is a
   * TIP quantity, d(dT_tip)/dt; the far-field reading is the tip value minus a
   * box-height constant, so any threshold on it is arbitrary rather than
   * physical.  Banking is better read post hoc from a windowed linear fit of
   * dT_tip against time.
   *
   * The far field is inert in this model (h'(0) = h''(0) = 0, noise strictly
   * interface-localised), so deep far-field undercooling is not itself a fault.
   *
   * The undercooling is still COMPUTED and LOGGED above; only the abort is gone.
   * The solid-cells-in-band abort above is retained: that one is a real failure
   * (the front has reached the top and the window can no longer track).
   */

  /* ---------------- short-box guards (quantities computed above) ------- */
  if (rank == MASTER) {
    printf("  short-box: tip row %ld (abs %ld, %ld cells below top), "
           "top-band dC/c0 = %+.3f%%, v_tip = ",
           tip_row, front_abs, (mpiparam.rows_y - 2) - tip_row,
           100.0 * dc_rel);
    if (have_vtip) {
      printf("%.4f mm/s = %.2f x V\n", vtip * 1.0e3,
             (V_pull != 0.0) ? vtip / V_pull : 0.0);
    } else {
      printf("(no previous sample)\n");
    }
    fflush(stdout);
  }

  if (dc_rel > dc_max_rel) {
    if (rank == MASTER) {
      printf("\nSHORT-BOX GUARD %s at step %ld: top-band composition is "
             "%+.3f%% vs c0 (limit %+.3f%%).\n"
             "  The solute boundary layer is no longer contained: with a "
             "no-flux top, rejected Si is being trapped and the alloy is\n"
             "  enriching.  The moving-window mass audit CANNOT see this -- its "
             "expected_after is defined to include refill at c0.\n"
             "  Increase the liquid above the tip (measure the decay distance; "
             "it is 22-40 l_D on a non-steady front, not 5-10).\n",
             shortbox_off ? "WARNING (abort disabled)" : "FAILED",
             absolute_step, 100.0 * dc_rel, 100.0 * dc_max_rel);
      fflush(stdout);
    }
    if (!shortbox_off) {
      MPI_Barrier(MPI_COMM_WORLD);
      MPI_Finalize();
      exit(SHORTBOX_EXIT_DC);
    }
  }

  if (tip_row >= 0 && ((mpiparam.rows_y - 2) - tip_row) < tip_margin_cells) {
    if (rank == MASTER) {
      printf("\nSHORT-BOX GUARD %s at step %ld: tip is %ld cells below the top "
             "boundary (limit %ld).\n"
             "  The front is running out of liquid.  Raise ny, or lower Shiftj "
             "so the moving window holds the front further down.\n",
             shortbox_off ? "WARNING (abort disabled)" : "FAILED",
             absolute_step, (mpiparam.rows_y - 2) - tip_row, tip_margin_cells);
      fflush(stdout);
    }
    if (!shortbox_off) {
      MPI_Barrier(MPI_COMM_WORLD);
      MPI_Finalize();
      exit(SHORTBOX_EXIT_TIP);
    }
  }

  if (have_vtip && V_pull > 0.0 && vtip > vtip_factor * V_pull) {
    if (rank == MASTER) {
      printf("\nSHORT-BOX GUARD %s at step %ld: v_tip = %.4f mm/s = %.2f x V "
             "(limit %.2f x).\n"
             "  This is the runaway signature: the front dumping "
             "banked undercooling into the liquid ahead of it.  It reached\n"
             "  15 x V there, and sum(phi) stayed at 1e-13 throughout.  Stop "
             "now; the frames after this point are not physics.\n",
             shortbox_off ? "WARNING (abort disabled)" : "FAILED",
             absolute_step, vtip * 1.0e3, vtip / V_pull, vtip_factor);
      fflush(stdout);
    }
    if (!shortbox_off) {
      MPI_Barrier(MPI_COMM_WORLD);
      MPI_Finalize();
      exit(SHORTBOX_EXIT_VTIP);
    }
  }
}

#endif
