
## Mesh commands

- [`POLY`](#poly)
- [`SEED`](#seed)
- [`QTRE`](#qtre)

### `POLY` mesh command <a name="poly"></a>

The `POLY` mesh command is used to define a heterogeneous solid, which is
described by its boundary and material interfaces. The boundary and material
interfaces are modeled as simple closed polygonal chains.


```
poly
num_poly, num_materials

// Definition: Ellipse
poly_nro, typ, material_nro, num_vertices, num_helpers, a, b, xc, yc, alpha

// Definition: Simple closed polygonal chain
poly_nro, typ, naterial_nro, num_vertices, num_helpers
    vertex_nro,num_divisions,x-coord,y-coord

    ...  // list vertices in counter clockwise order

    num_vertices

    vertex_nro,x-coord,y-coord

    ...

    num_helpers

...

num_poly
```

| Parameter | Description |
| :--- | :--- |
| num_poly | number of polygons |
| num_materials | number of material sets|
| poly_nro | polygon number |
| typ | polygon typ |
| material_nro | material set number |
| num_vertices | number of vertices |
| num_helpers | number of helperpoints |
| a | semi-major axis |
| b | semi-minor axis |
| xc | x-coord ellipse center |
| yc | y-coord ellipse xenter |
| alpha | rotate ellipse $\alpha$ degrees counterclockwise with respect to the x-axis |
| vertex_nro|  vertex number |
|num_divisions| number of divisions |
|x-coord| vertex/helperpoint x-coord |
|y-coord| vertex/helperpoint y-coord |


### `SEED` mesh command <a name="seed"></a>
The `SEED` mesh command is used to define seeding points. The seeding points
can be used to control the local density of the quadtree mesh.

```
seed (optional)
num_seeds
seed_nro,num_divisions,x-ccord,y-coord

...

num_seeds
```



| Parameter | Description |
| :--- | :--- |
| num_seeds | number of seeding points |
| seed_nro |  seeding point number |
| num_divisions | number of divisions |
| x-coord | seeding point x-coord |
| y-coord | seeding point y-coord |

### `QTRE` mesh commad <a name="qtre"></a>
The `QTRE` mesh command is used to specify the parameters to control the quadtree
decomposition


```
qtree
level_min, max_seeds_cell

```


| Parameter | Description |
| :--- | :--- |
| level_min |  minimum division level |
| max_seeds_cell |  maximum seeding points per cell |



