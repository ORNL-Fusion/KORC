
# KORC


 the Kinetic Orbit Runaway electrons Code (KORC) follows relativistic electrons in general electric and magnetic fields under the full Lorentz force, collisions, and radiation losses.


## Quick-Start

1. Find your platform under the `.cmfkit/platform` directory

```bash
[user@localhost korc]$ tree .cmfkit/platform/

.cmfkit/platform/
├── github-vm-macos-latest
│   ├── match
│   └── setup
├── github-vm-ubuntu-latest
│   ├── match
│   └── setup
├── github-vm-windows-latest
│   ├── match
│   └── setup
├── metal-macos-12
│   ├── match
│   └── setup
├── readybake
│   └── Dockerfile
├── readybake-alpine
│   └── Dockerfile
├── readybake-centos
│   └── Dockerfile
├── readybake-debian
│   └── Dockerfile
└── readybake-ubuntu
    └── Dockerfile

```


2. Create a `cmfkit_vars` file at the base directory of the project and set the `CMFKIT_PLATFORM` variable to match the platform name.

```bash

CMFKIT_PLATFORM='metal-macos-12'
CMFKIT_SCENARIO='basic'
```


3. Run the launcher script under the .cmfkit_base directory

```
[user@localhost korc]$ ./.cmfkit_base/launch build
```


## More Info

http://ornl-fusion.github.io/KORC

