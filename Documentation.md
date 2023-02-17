
## Mesh commands

### POLY mesh command

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


|   |   |
| --- | --- |
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
| alpha | rotate ellipse alpha degrees counterclockwise with respect to the x-axis |  
| vertex_nro|  vertex number |
|num_divisions| number of divisions |
|x-coord| vertex/helperpoint x-coord |
|y-coord| vertex/helperpoint y-coord |


### SEED mesh command

### QTREE mesh commad
