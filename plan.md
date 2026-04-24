right, so we have solved the issue of sherds not overlapping each other in the sense that they can't be phyiscally inside each other as that will trigger the points penalty that are colliding.

The only problems remains now is that when two sherds occupy the same position but slightly on top of each other so that they don't trigger the collision/overlap penalty.

A rough idea for detecting this penalty would be to get the rotation and position of both sherds and compare them to see if they are similar. if yes, then apply the penalty.

what do you think? don't code just think and let me know your thoughts.



Your idea is good, but I would not use pose similarity alone.

1. Why your idea is valid
1. If two sherds are almost in the same pose, relative transform is close to identity:
$$
T_{ij}=T_j^{-1}T_i,\quad
t_{ij}=\|\Delta t_{ij}\|,\quad
\theta_{ij}=\text{angle}(\Delta R_{ij})
$$
2. Small $t_{ij}$ and small $\theta_{ij}$ are exactly the signature of stacking/duplicate placement.

2. Why pose-only can be risky
1. Two correctly placed neighboring sherds can also have similar orientation.
2. So a pure pose penalty may punish valid local assembly.

3. Best practical design
Use a gated near-coincidence penalty:
1. Pose gate:
$$
g_{pose}=\max\!\left(0,1-\frac{t_{ij}}{\tau_t}\right)\cdot
\max\!\left(0,1-\frac{\theta_{ij}}{\tau_r}\right)
$$
2. Near-band cloud coverage (not overlap): fraction of points with nearest-neighbor distance in
$$
\epsilon < d_{nn} < \epsilon_{near}
$$
3. Final added term:
$$
P_{stack}=\lambda_{stack}\,g_{pose}\,\phi_{near}
$$
where $\phi_{near}$ is near-band coverage ratio.

This catches “almost on top but not penetrating” and avoids penalizing tiny boundary contacts.

4. Why this is better than overlap-only in your code
Your current overlap term in [class/genetic_algorithm.h](class/genetic_algorithm.h#L1708) only triggers for $d<\epsilon$, and returns zero when there are no hits in genetic_algorithm.h.  
So stacked-but-just-outside-epsilon can slip through.

5. Tuning guidance
1. Set $\tau_t$ relative to sherd size, not fixed global units.
2. Set $\tau_r$ modest (around small-angle similarity only).
3. Keep $\epsilon_{near}$ larger than collision epsilon and also meaningful relative to voxel size.
4. Start with small $\lambda_{stack}$, then raise until these failure cases disappear.

So overall: yes, add pose similarity, but only as a gate multiplied by a near-band coverage metric. That will target exactly the failure mode you described.