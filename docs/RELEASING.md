# Releasing AlgRemez.jl

## Before the first registration

1. Publish the repository under the name `AlgRemez.jl`.
2. Keep the package name `AlgRemez` and UUID
   `12fcd701-c184-4ee5-aece-092528ad8aba` unchanged.
3. Confirm that CI passes on Julia 1.10 and the latest stable Julia release.
4. Confirm that the repository is public and that `LICENSE` is present.
5. Install the Julia Registrator GitHub App or use its web interface.
6. Trigger the initial registration from a commit on the default branch.

The package source must not include the historical C++ implementation.
`AlgRemez_jll` remains a test-only dependency and is not loaded by package
users.

## Version releases

1. Update `version` in `Project.toml` according to semantic versioning.
2. Add the release notes to `CHANGELOG.md`.
3. Run the full test suite and the RHMC benchmark.
4. Commit and push the release preparation.
5. Trigger Registrator for that commit.
6. Let TagBot create the matching Git tag and GitHub release after registry
   acceptance.

Do not reuse a version number after changing a registered source tree; increase
the patch version instead.
