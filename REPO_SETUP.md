# How to put this on GitHub

The git repository is already initialised here with one commit on the
`main` branch. You only need to create the empty repository on GitHub and
push.

## 1. Create the repository

On github.com: New repository.

- Name: `fluent-udf-library`
- Description (paste this):

      Tested User Defined Functions for ANSYS Fluent: rigid-body and six-DOF motion, VOF numerical wave tank, boundary profiles and material properties.

- Public
- Do NOT tick "Add a README", "Add .gitignore" or "Choose a license".
  All three are already in the commit, and ticking them creates a
  conflicting first commit that you would then have to merge.

## 2. Push

Copy the URL GitHub shows you, then from the folder that contains this
file:

    git remote add origin https://github.com/YOUR_USERNAME/fluent-udf-library.git
    git push -u origin main

If it asks for a password, GitHub wants a personal access token, not your
account password: Settings > Developer settings > Personal access tokens >
Fine-grained tokens, with Contents: Read and write on this repository.

## 3. Topics

Settings, or the gear next to "About" on the repository page. Topics are
what makes it findable in GitHub search:

    ansys, ansys-fluent, fluent, udf, cfd, computational-fluid-dynamics,
    wave-tank, vof, six-dof, dynamic-mesh, marine-engineering,
    free-surface, hydrodynamics, c

## 4. Two settings worth changing

- About panel: tick "Releases" and "Packages" off if you find them noisy;
  leave Issues on.
- Settings > General > Features: make sure Issues is enabled, otherwise
  the templates in `.github/ISSUE_TEMPLATE/` do nothing.

## 5. A release (optional, later)

When you have run a couple of these in Fluent and are happy, tag it:

    git tag -a v0.1.0 -m "First release"
    git push origin v0.1.0

Then Releases > Draft a new release > pick the tag. A tagged release is
what people cite and what Zenodo picks up if you ever want a DOI for it.

## Working after that

    git add -A
    git commit -m "what changed"
    git push

Run `sh tests/run_tests.sh` before each push. It takes a few seconds and
catches the kind of mistake that is embarrassing in a public repository.
