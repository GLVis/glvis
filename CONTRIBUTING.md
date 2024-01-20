<p align="center">
<a href="https://glvis.org/"><img alt="glvis" src="https://glvis.org/img/logo-300.png"></a>
</p>

<p align="center">
<a href="https://github.com/GLVis/glvis/releases/latest"><img alt="Release" src="https://img.shields.io/badge/release-v4.2-success.svg"></a>
<a href="https://github.com/GLVis/glvis/actions/workflows/builds.yml"><img alt="Build" src="https://github.com/GLVis/glvis/actions/workflows/builds.yml/badge.svg"></a>
<a href="https://github.com/glvis/glvis/blob/master/LICENSE"><img alt="License" src="https://img.shields.io/badge/License-BSD-success.svg"></a>
<a href="https://glvis.github.io/doxygen/html/index.html"><img alt="Doxygen" src="https://img.shields.io/badge/code-documented-success.svg"></a>
<a href="https://glvis.github.io/releases/glvis-macOS.dmg"><img alt="License" src="https://img.shields.io/badge/Download-Mac-success.svg"></a>
<a href="https://glvis.github.io/releases/glvis-Windows.zip"><img alt="License" src="https://img.shields.io/badge/Download-Windows-success.svg"></a>
</p>


# How to Contribute

The GLVis team welcomes contributions at all levels: bugfixes; code improvements;
simplifications; new visualization capabilities; improved documentation; etc.

GLVis is distributed under the terms of the BSD-3 license. All new contributions
must be made under this license.

If you plan on contributing to GLVis, consider reviewing the
[issue tracker](https://github.com/glvis/glvis/issues) first to check if a thread
already exists for your desired feature or the bug you ran into. Use a pull
request (PR) toward the `glvis:master` branch to propose your contribution. If
you are planning significant code changes or have questions, you may want to
open an [issue](https://github.com/glvis/glvis/issues) before issuing a PR. In
addition to technical contributions, we are also interested in your results and
[simulation images](https://glvis.org/gallery/), which you can share via a pull
request in the [glvis/web](https://github.com/glvis/web) repo.

See the [Quick Summary](#quick-summary) section for the main highlights of our
GitHub workflow. For more details, consult the following sections and refer
back to them before issuing pull requests:

- [Quick Summary](#quick-summary)
- [Code Overview](#code-overview)
- [GitHub Workflow](#github-workflow)
  - [GLVis Organization](#glvis-organization)
  - [New Feature Development](#new-feature-development)
  - [Developer Guidelines](#developer-guidelines)
  - [Pull Requests](#pull-requests)
  - [Pull Request Checklist](#pull-request-checklist)
  - [Releases](#releases)
  - [Release Checklist](#release-checklist)
- [LLNL Workflow](#llnl-workflow)
- [Automated Testing](#automated-testing)
- [Contact Information](#contact-information)

Contributing to GLVis requires knowledge of Git and, likely, OpenGL and/or finite
elements. If you are new to Git, see the [GitHub learning
resources](https://help.github.com/articles/git-and-github-learning-resources/).
To learn more about the finite element method, see this [MFEM page](https://mfem.org/fem).

*By submitting a pull request, you are affirming the [Developer's Certificate of
Origin](#developers-certificate-of-origin-11) at the end of this file.*


## Quick Summary

- We encourage you to [join the GLVis organization](#glvis-organization) and create
  development branches off `glvis:master`.
- Please follow the [developer guidelines](#developer-guidelines), in particular
  with regards to documentation and code styling.
- Pull requests  should be issued toward `glvis:master`. Make sure
  to check the items off the [Pull Request Checklist](#pull-request-checklist).
- When your contribution is fully working and ready to be reviewed, add
  the `ready-for-review` label.
- PRs are treated similarly to journal submission with an "editor" assigning two
  reviewers to evaluate the changes.
- The reviewers have 3 weeks to evaluate the PR and work with the author to
  fix issues and implement improvements.
- We use [milestones](https://github.com/glvis/glvis/milestones) to coordinate the
  work on different PRs toward a release.
- Don't hesitate to [contact us](#contact-information) if you have any questions.


### Code Overview

#### Source code structure:

The GLVis library uses object-oriented design principles which reflect, in code,
the independent mathematical concepts of meshing, linear algebra and finite
element spaces and operators.

The GLVis source code has the following structure:

```
  .
  ├── lib
  │   └── gl
  │       └── shaders
  ├── share
  └── tests
```

## GitHub Workflow

The GLVis GitHub organization: https://github.com/glvis, is the main developer hub
for the GLVis project.

If you plan to make contributions or would like to stay up-to-date with changes
in the code, *we strongly encourage you to [join the GLVis organization](#glvis-organization)*.

This will simplify the workflow (by providing you additional permissions), and
will allow us to reach you directly with project announcements.


### GLVis Organization

#### Getting started (GitHub)
Before you can start, you need a GitHub account, here are a few suggestions:
  + Create the account at: [github.com/join](https://github.com/join).
  + For easy identification, please add your full name and maybe a picture of you at:
    https://github.com/settings/profile.
  + To receive notification, set a primary email at: https://github.com/settings/emails.
  + For password-less pull/push over SSH, add your SSH keys at: https://github.com/settings/keys.

#### Joining
- [Contact us](#contact-information) for an invitation to join the GLVis GitHub
  organization. You will receive an invitation email, which you can directly accept.
  Alternatively, *after logging into GitHub*, you can accept the invitation at
  the top of https://github.com/glvis.
- Consider making your membership public by going to https://github.com/orgs/glvis/people
  and clicking on the organization visibility drop box next to your name.
- Project discussions and announcements will be posted at
  https://github.com/orgs/glvis/teams/everyone.

#### Structure
- The GLVis source code is in the [glvis](https://github.com/glvis/glvis)
  repository.
- The website and corresponding documentation are in the
  [web](https://github.com/glvis/web) repository.



### New Feature Development

- A new feature should be important enough that at least one person, the
  author, is willing to work on it and be its champion.

- The author creates a branch for the new feature (with suffix `-dev`), off
  the `master` branch, or another existing feature branch, for example:

  ```
  # Clone assuming you have setup your ssh keys on GitHub:
  git clone git@github.com:glvis/glvis.git

  # Alternatively, clone using the "https" protocol:
  git clone https://github.com/glvis/glvis.git

  # Create a new feature branch starting from "master":
  git checkout master
  git pull
  git checkout -b feature-dev

  # Work on "feature-dev", add local commits
  # ...

  # (One time only) push the branch to github and setup your local
  # branch to track the github branch (for "git pull"):
  git push -u origin feature-dev

  ```

- **We prefer that you create the new feature branch inside the GLVis organization
  as opposed to in a fork.** This allows everyone in the community to collaborate
  in one central place.

  - If you prefer to work in your fork, please [enable upstream edits](https://help.github.com/articles/allowing-changes-to-a-pull-request-branch-created-from-a-fork/).

  - Never use the `next` branch to start a new feature branch!

- The typical feature branch name is `new-feature-dev`, e.g. `pumi-dev`. While
  not frequent in GLVis, other suffixes are possible, e.g. `-fix`, `-doc`, etc.


### Developer Guidelines

- *Keep the code lean and as simple as possible*
  - Well-designed simple code is frequently more general and powerful.
  - Lean code base is easier to understand by new collaborators.
  - New features should be added only if they are necessary or generally useful.
  - Introduction of language constructions not currently used in GLVis should be
    justified and generally avoided (to maintain portability to various systems
    and compilers, including early access hardware).
  - We prefer basic C++ and the C++03 standard, to keep the code readable by
    a large audience and to make sure it compiles anywhere.

- *Keep the code general and reasonably efficient*
  - The main goal is fast prototyping for research and application development.
  - When in doubt, generality wins over efficiency.
  - Respect the needs of different users (current and/or future).

- *Keep things separate and logically organized*
  - General usage features go in GLVis (implemented in as much generality as
    possible), non-general features go into external apps.
  - Inside GLVis, compartmentalize between linalg, fem, mesh, GLVis, etc.
  - Contributions that are project-specific or have external dependencies are
    allowed (if they are of broader interest), but should be `#ifdef`-ed and not
    change the code by default.

- Code specifics
  - All significant new classes, methods and functions have Doxygen-style
    documentation in source comments.
  - Consistent code styling is enforced with `make style` in the top-level
    directory. This requires [Artistic Style](http://astyle.sourceforge.net)
    version 3.1 and MFEM's style configuration file, typically located in
    `../mfem/config/mfem.astylerc`.
  - When manually resolving conflicts during a merge, make sure to mention the
    conflicted files in the commit message.


### Pull Requests

- When your branch is ready for other developers to review / comment on
  the code, create a pull request towards `glvis:master`.

- Pull request typically have titles like:

     `Description [new-feature-dev]`

  for example:

     `Parallel Unstructured Mesh Infrastructure (PUMI) integration [pumi-dev]`

  Note the branch name suffix (in square brackets).

- Titles may contain a prefix in square brackets to emphasize the type of PR.
  Common choices are: `[DON'T MERGE]`, `[WIP]` and `[DISCUSS]`, for example:

     `[DISCUSS] Hybridized DG [hdg-dev]`

- If the PR is still a work in progress, add the `WIP` label to it and
  optionally the `[WIP]` prefix in the title.

- Add a description, appropriate labels and assign yourself to the PR. The GLVis
  team will add reviewers as appropriate.

- List outstanding TODO items in the description, see PR #222 for an example.

- When your contribution is fully working and ready to be reviewed, add
  the `ready-for-review` label.

- PRs are treated similarly to journal submission with an "editor" assigning
  two reviewers to evaluate the changes. The reviewers have 3 weeks to evaluate
  the PR and work with the author to implement improvements and fix issues.

### Pull Request Checklist

Before a PR can be merged, it should satisfy the following:

- [ ] Code builds.
- [ ] Code passes `make style`.
- [ ] Update `CHANGELOG`:
    - [ ] Is this a new feature users need to be aware of?
    - [ ] Does it make sense to create a new section in the `CHANGELOG` to group with other related features?
- [ ] Update `INSTALL`:
    - [ ] Had a new optional library been added? If so, what range of versions of this library are required? (*Make sure the external library is compatible with our BSD license, e.g. it is not licensed under GPL!*)
    - [ ] Have the version ranges for any required or optional libraries changed?
    - [ ] Does `make` or `cmake` have a new target?
    - [ ] Did the requirements or the installation process change? *(rare)*
- [ ] Update `.gitignore`:
    - [ ] Check if `make distclean; git status` shows any files that were generated from the source by the project (not an IDE) but we don't want to track in the repository.
    - [ ] Add new patterns (just for the new files above) and re-run the above test.
- [ ] New capability:
   - [ ] All significant new classes, methods and functions have Doxygen-style documentation in source comments.
   - [ ] Consider saving cool simulation pictures with the new capability in the Confluence gallery (LLNL only) or submitting them, via pull request, to the gallery section of the `glvis/web` repo.
   - [ ] If this is a major new feature, consider mentioning it in the short summary inside `README` *(rare)*.
   - [ ] List major new classes in `doc/CodeDocumentation.dox` *(rare)*.
- [ ] Update this checklist, if the new pull request affects it.



### Releases

- Releases are just tags in the `master` branch, e.g. https://github.com/glvis/glvis/releases/tag/v3.4,
  and have a version that ends in an even "patch" number, e.g. `v3.4.2` or
  `v3.4` (by convention `v3.4` is the same as `v3.4.0`.)  Between releases, the
  version ends in an odd "patch" number, e.g. `v3.4.1`.

- We use [milestones](https://github.com/glvis/glvis/milestones) to coordinate the
  work on different PRs toward a release, see for example the
  [v3.4 release](https://github.com/glvis/glvis/milestone/1?closed=1).

- After a release is complete, the `next` branch is recreated, e.g. as follows
  (replace `3.4` with current release):
  - Rename the current `next` branch to `next-pre-v3.4`.
  - Create a new `next` branch starting from the `v3.4` release.
  - Local copies of `next` can then be updated with `git fetch origin next && git checkout -B next origin/next`.

### Release Checklist
- [ ] Update the GLVis version in the following files:
    - [ ] `CHANGELOG`
    - [ ] `README.md`
    - [ ] `vcpkg.json`
    - [ ] `share/Info.plist`
    - [ ] `share/Info.cmake.plist.in`
- [ ] Check that version requirements for each of GLVis's dependencies are documented in `INSTALL` and up-to-date
- [ ] Update the `CHANGELOG` to organize all release contributions
- [ ] Review the whole source code once over
- [ ] Ask GLVis-based applications to test the pre-release branch
- [ ] Test on additional platforms and compilers
- [ ] Tag the repository:

  ```
  git tag -a v3.1 -m "Official release v3.1"
  git push origin v3.1
  ```
- [ ] Create the release tarball and push to `glvis/releases`.
- [ ] Recreate the `next` branch as described in previous section.
- [ ] Update and push documentation  to `glvis/doxygen`.
- [ ] Update URL shortlinks:
    - [ ] Create a shortlink at [http://bit.ly/](http://bit.ly/) for the release tarball, e.g. https://glvis.github.io/releases/glvis-3.1.tgz.
    - [ ] (LLNL only) Add and commit the new shortlink in the `links` and `links-glvis` files of the internal `glvis/downloads` repo.
    - [ ] Add the new shortlinks to the GLVis packages in `spack`, `homebrew/science`, `VisIt`, etc.
- [ ] Update website in `glvis/web` repo:
    - Update version and shortlinks in `src/index.md` and `src/download.md`.
    -  Use [cloc-1.62.pl](http://cloc.sourceforge.net/) and `ls -lh` to estimate the SLOC and the tarball size in `src/download.md`.


## LLNL Workflow

### Mirroring on Bitbucket

- The GitHub `master` and `next` branches are mirrored to the LLNL institutional
  Bitbucket repository as `gh-master` and `gh-next`.

- `gh-master` is merged into LLNL's internal `master` through pull requests; write
  permissions to `master` are restricted to ensure this is the only way in which it
  gets updated.

- We never push directly from LLNL to GitHub.

- Versions of the code on LLNL's internal server, from most to least stable:
  - GLVis official release on glvis.org -- Most stable, tested in many apps.
  - `glvis:master` -- Recent development version, guaranteed to work.
  - `glvis:gh-master` -- Stable development version, passed testing, you can use
     it to build your code between releases.
  - `glvis:gh-next` -- Bleeding-edge development version, may be broken, use at
     your own risk.


## Contact Information

- Contact the GLVis team by posting to the [GitHub issue tracker](https://github.com/glvis/glvis).
  Please perform a search to make sure your question has not been answered already.



## [Developer's Certificate of Origin 1.1](https://developercertificate.org/)

By making a contribution to this project, I certify that:

(a) The contribution was created in whole or in part by me and I have the right
    to submit it under the open source license indicated in the file; or

(b) The contribution is based upon previous work that, to the best of my
    knowledge, is covered under an appropriate open source license and I have
    the right under that license to submit that work with modifications, whether
    created in whole or in part by me, under the same open source license
    (unless I am permitted to submit under a different license), as indicated in
    the file; or

(c) The contribution was provided directly to me by some other person who
    certified (a), (b) or (c) and I have not modified it.

(d) I understand and agree that this project and the contribution are public and
    that a record of the contribution (including all personal information I
    submit with it, including my sign-off) is maintained indefinitely and may be
    redistributed consistent with this project or the open source license(s)
    involved.
