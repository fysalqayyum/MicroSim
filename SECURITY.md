# Security and sensitive-data policy

Do not open a public issue containing a credential or private dataset. Report a
suspected repository secret privately to the fork maintainer through the
contact method on the maintainer's GitHub profile.

The repository must not contain:

- passwords, API tokens, private keys, certificates, or credential files;
- personal machine paths, usernames, SSH aliases, or scheduler accounts added
  by this fork;
- private thermodynamic databases or data with unclear redistribution terms;
- raw production outputs that may carry confidential project metadata.

If a secret is committed, deleting it in a later commit is insufficient.
Revoke or rotate it first, then remove it from the complete Git history before
publishing the corrected branch.

Scientific defects, numerical instability, or incorrect documentation can be
reported through ordinary GitHub issues when no sensitive data are attached.
