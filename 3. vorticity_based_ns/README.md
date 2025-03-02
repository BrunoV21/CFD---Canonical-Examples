## Vorticity Based Navier Stokes 

[](./assets/u_magnitude.png)

Based on Joe Molvar's work and A. Salih notes, I implemented  a simple vorticity based NS solver. The advantage of the vorticity stream function formulation is that we do not need to solve for pressure so we can geet the 2D flow field by just solving two equations.

You can refer to thepds available in the resources dir.

Here are a couple of mages of flow over a backward facing step.

[](./assets/streamlines.mp4)

[](./assets/code_verification.png)