from typing import List, Tuple, Sequence, Iterable
from jinja2 import Template
from IPython.core.display import HTML
from numpy import ndarray
import numpy
from io import StringIO
import time
from collections import OrderedDict
import os

# local paths to javascript and css files
CURRENT_PATH = os.path.dirname(os.path.realpath(__file__))
STATIC_ROOT = os.path.join(CURRENT_PATH, 'static')
X3D_ROOT = os.path.join(STATIC_ROOT, 'x3d')
# X3DOM_JS = os.path.join(X3D_ROOT, 'x3dom.debug.js')  # 'x3dom.js')
# X3DOM_CSS = os.path.join(X3D_ROOT, 'x3dom.css')

X3DOM_JS = 'https://x3dom.org/download/dev/x3dom.js'
X3DOM_CSS = 'https://x3dom.org/download/dev/x3dom.css'

COLORS_RGB = OrderedDict([
    ('pure magenta red', ((237,  20,  91), '#ed145b')),
    ('pure magenta', ((236, 0, 140), '#ec008c')),
    ('pure violet magenta', ((146, 39, 143), '#92278f')),
    ('pure violet', ((102, 45, 145), '#662d91')),
    ('pure blue violet', ((46, 49, 146), '#2e3192')),
    ('pure blue', ((0, 84, 166), '#0054a6')),
    ('pure cyan blue', ((0, 114, 188), '#0072bc')),
    ('pure cyan', ((0, 174, 239), '#00aeef')),
    ('pure green cyan', ((0, 168, 157), '#00a99d')),
    ('pure green', ((0, 166, 81), '#00a651')),
    ('pure yellow green', ((57, 181, 74), '#39b54a')),
    ('pure pea green', ((141, 198, 63), '#8dc63f')),
    ('pure yellow', ((255, 242, 0), '#fff200')),
    ('pure yellow orange', ((247, 148, 29), '#f7941d')),
    ('pure red orange', ((242, 101, 34), '#f26522')),
    ('pure red', ((237, 28, 36), '#ed1c24'))
])


def bounding_box(coordinate_arrays: Sequence[ndarray]) -> ndarray:
    """
    A bounding box containing all the coordinates in a collection of coordinate arrays.
    
    Args:
        coordinate_arrays: collection of :math:`i = 1, ..., N` arrays, each with shape :math:`(M_i, N)`, or :math:`M_i`
            points, where :math:`N` is the coordinate dimension, e.g., 2 or 3 for 2D or 3D coordinates.

    Returns:
        bbox - array containing minimum and maximum coordinates with shape (2, NDIM) 
    """
    # get a bounding box for each of the arrays
    a_mins = list()
    a_maxs = list()
    for a in coordinate_arrays:
        a_mins.append(numpy.amin(a, axis=0))
        a_maxs.append(numpy.amax(a, axis=0))
    # find the bounding box containing all the bounding boxes
    b_min = a_mins[0]
    b_max = a_maxs[0]
    for i in range(1, len(coordinate_arrays)):
        numpy.minimum(b_min, a_mins[i], b_min)
        numpy.maximum(b_max, a_maxs[i], b_max)
    return numpy.array([b_min, b_max])


def normalize_v3(vectors: ndarray) -> ndarray:
    """
    Normalize a numpy array of 3 component vectors in place.

    Args:
        vectors: the input array of shape=(N,3)

    Returns:
        vectors - normalized vectors array
    """
    lens = numpy.sqrt(vectors[:, 0] ** 2 + vectors[:, 1] ** 2 + vectors[:, 2] ** 2)
    vectors[:, 0] /= lens
    vectors[:, 1] /= lens
    vectors[:, 2] /= lens
    return vectors


def triangle_normals(indices: ndarray, points: ndarray, normal_per_vertex=False) -> ndarray:
    """
    Compute the normals for triangles given by indices into an array of vertex coordinates.

    The normal is the vector of length 1.0 using the direction defined by the CCW vertices.

    Args:
        indices: triangle vertices given by index into the points array, with shape(N, 3)
        points: array containing vertex coordinates with shape(M, 3)
        normal_per_vertex: if True, average the normals from each triangle sharing a vertex

    Returns:
        normals - array of one normal vector per triangle with shape(N,3), or if normal_per_vertex is True
            an array of one normal vector per vertex with shape(M,3), with dtype=points.dtype

    If a vertex coordinate is not referenced in the indices array and normal_per_vertex is True,
    the corresponding entry in the normals array will be the zero vector.

    Raises:
        IndexError: if any entries in indices >= M
    """

    # create an indexed view into the points array using the array of three indices for triangles
    tris = points[indices]

    # vector for each triangle, by taking the cross product of the vectors p1-p0, and p2-p0
    n = numpy.cross(tris[::, 1] - tris[::, 0], tris[::, 2] - tris[::, 0])
    normalize_v3(n)

    # now we have a normalized array of normals, one per triangle (for flat shading).
    # are we going to average to vertices?
    if normal_per_vertex:
        # Multiple triangles contribute to every vertex; add them up and normalize again.
        normals = numpy.zeros(points.shape, dtype=points.dtype)
        num_faces = indices.shape[0]
        for local_vertex in range(3):
            for global_face in range(num_faces):
                normals[indices[global_face, local_vertex]] += n[global_face]
        normalize_v3(normals)
        return normals
    else:
        return n


def notebook_scene(name_space: str, head: str, subhead: str, description: str, width: int, height: int, scene_file: str,
                   show_log: bool = False, center: Iterable = (0.0, 0.0, 0.0)) -> HTML:
    """
    Display a 3D model inside an ipython notebook.

    Return valid HTML code that will display interactive 3D model using the X3DOM framework. The the model is
    given as a x3d scene using the schema at `http://www.web3d.org/specifications/x3d-3.3.xsd`.

    Example:

       Inside a Jupyter notebook cell:

       >>> notebook_scene('Model_0', 'Model', 'A solid mechanics model', 'series 3c', 500, 500, 'myModel.x3d')

    Args:
        name_space: unique short name identifier for avoiding conflicts with multiple scenes
        head: title for the model
        subhead: subtitle for the model
        description: longer text description for this model view
        width: width of the display window in pixels
        height: height of the display window in pixels
        scene_file: an existing scene file name
        show_log: if True, log template expansion
        center: center of rotation

    Returns:
        html code block for embedding inside a <body> tag or as part of a ipython Jupyter notebook.
    """
    if show_log:
        show_log_str = 'true'
    else:
        show_log_str = 'false'

    t = Template(WRAPPER_TEMPLATE)
    source = t.render(x3dom_js=X3DOM_JS, x3dom_css=X3DOM_CSS, namespace=name_space, head=head, subhead=subhead,
                      description=description, show_log=show_log_str, width=width, height=height,
                      scene_file=scene_file, center=','.join(map(str, center)))
    return HTML(source)


def notebook_gltf_box() -> HTML:
        """
        Display a gltf 3D model box in an ipython notebook.

        Return valid HTML code that will display interactive 3D model using the X3DOM framework.

        Returns:
            html code block for embedding inside a <body> tag or as part of a ipython Jupyter notebook.
        """

        t = Template(X3D_GLTF_TEMPLATE)
        source = t.render(x3dom_js=X3DOM_JS)
        return HTML(source)


def write_scene(file: str, title: str, summary: str, faces: List[Tuple[ndarray, ndarray]], use_binary: bool = False):
    """
    Write X3D scene file.

    Uses the schema at `http://www.web3d.org/specifications/x3d-3.1.xsd`.

    Args:
        use_binary: use a binary format file instead of plain text
        file: path to file name, an ".x3d" extension will be appended to it
        title: title of scene
        summary: description of origin
        faces: list of pairs of face to vertex connectivity and vertex coordinates (face_to_vertex, coordinates)

    """
    if use_binary:
        directory = os.path.dirname(os.path.abspath(file))
        shapes_buffer = _create_shapes_buffer_binary(directory, faces)
    else:
        shapes_buffer = _create_shapes_buffer(faces)
    try:
        func_name = 'affect.display.write_scene'
        now = time.strftime("%H:%M:%S")
        date = time.strftime("%Y-%m-%d")
        t = Template(MODEL_TEMPLATE)
        result = t.render(file=file,
                          title=title,
                          function=func_name,
                          summary=summary,
                          date=date,
                          time=now,
                          shapes=shapes_buffer.getvalue())
    except Exception:
        raise RuntimeError('Could not render scene to x3d file')
    else:
        shapes_buffer.close()  # make certain we close the buffer
    with open(file + '.x3d', 'w') as text_file:
        text_file.write(result)


def _create_shapes_buffer_binary(directory: str, face_blocks: List[Tuple[ndarray, ndarray]]) -> StringIO:
    # noinspection PyPep8
    """
    Write multiple face blocks and coordinates to the Shape elements of an x3d scene.

    Args:
        face_blocks: list of pairs of face to vertex connectivity and vertex coordinates (face_to_vertex, coordinates)

    Returns:
        string buffer of IndexedFaceSet items of the form (this one is a hexahedron):

        <IndexedFaceSet coordIndex='0 1 5 4 -1 1 2 6 5 -1 2 3 7 6 -1 3 0 4 7 -1 0 3 2 1 -1 4 5 6 7 -1'>
          <Coordinate point='-1.0 -1.0 1.0 1.0 -1.0 1.0 1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 -1.0 -1.0 1.0 -1.0'/>
        </IndexedFaceSet>
    """
    shapes_buffer = StringIO()
    shape_t = Template(SHAPE_TEMPLATE)
    binary_triangles_t = Template(BINARY_GEOMETRY_TEMPLATE)

    color = "0.8 0.8 0.95"
    shape_count = 0

    index_name = 'index.bin'
    index_path = os.path.join(directory, index_name)
    coord_name = 'coord.bin'
    coord_path = os.path.join(directory, coord_name)
    normal_name = 'normal.bin'
    normal_path = os.path.join(directory, normal_name)

    for face_to_vertex, coordinates in face_blocks:

        with open(index_name, 'wb') as index_file:
            # face_to_vertex.astype(numpy.uint16).tofile(index_path)
            # vertex_count = face_to_vertex.size

            indices = numpy.array([[0, 1, 2], [1, 3, 2], [0, 2, 3], [0, 3, 1]])
            indices.astype(numpy.uint16).tofile(index_path)
            vertex_count = indices.size

        with open(coord_name, 'wb') as coord_file:

            # coordinates.astype(numpy.float32).tofile(coord_path)

            coords = numpy.array([[0.5, 0.0, 1.0], [1.0, 0.0, 0.5], [0.5, 1.0, 0.5], [0.0, 0.0, 0.0]])
            coords.astype(numpy.float32).tofile(coord_path)

        # normals = average_normals_per_vertex(face_to_vertex, coordinates)
        normals = triangle_normals(indices, coords)
        print('******* normals =\n{}'.format(normals))

        with open(normal_name, 'wb') as normal_file:
            normals.astype(numpy.float32).tofile(normal_path)

        binary_geometry = binary_triangles_t.render(vertex_count=vertex_count,
                                                    index=index_name,
                                                    coord=coord_name,
                                                    normal=normal_name)

        shape_result = shape_t.render(shape_name=str(shape_count),
                                      color=color,
                                      geometry_nodes=binary_geometry)

        shapes_buffer.write(shape_result)
        shape_count = shape_count + 1

        break   # ********* for experiment with binary file, can't get flat face normals to display properly

    return shapes_buffer


def _create_shapes_buffer(face_blocks: List[Tuple[ndarray, ndarray]]) -> StringIO:
    # noinspection PyPep8
    """
    Write multiple face blocks and coordinates to the Shape elements of an x3d scene in textual format.

    This version of _create_shapes_buffer_binary is better suited for debugging since it writes data inline in the
    x3dom file in text format as opposed to binary array data.

    Args:
        face_blocks: list of pairs of face to vertex connectivity and vertex coordinates (face_to_vertex, coordinates)

    Returns:
        string buffer of IndexedFaceSet items of the form (this one is a hexahedron):

        <IndexedFaceSet coordIndex='0 1 5 4 -1 1 2 6 5 -1 2 3 7 6 -1 3 0 4 7 -1 0 3 2 1 -1 4 5 6 7 -1'>
          <Coordinate point='-1.0 -1.0 1.0 1.0 -1.0 1.0 1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 -1.0 -1.0 1.0 -1.0'/>
        </IndexedFaceSet>
    """
    shapes_buffer = StringIO()
    shape_t = Template(SHAPE_TEMPLATE)
    triangle_set_t = Template(INDEXED_TRIANGLE_SET_TEMPLATE)

    color = "0.8 0.8 0.95"
    shape_count = 0

    for face_to_vertex, coordinates in face_blocks:
        # write the unstructured face to vertex, and vertex coordinates
        coord_index = ' '.join(map(str, face_to_vertex.flat))
        point = ' '.join(map(lambda x: '{:.6e}'.format(x) if x > 100 else '{:.6f}'.format(x), coordinates.flat))
        shape_result = shape_t.render(shape_name=str(shape_count),
                                      color=color,
                                      geometry_nodes=triangle_set_t.render(coord_index=coord_index, point=point))
        shapes_buffer.write(shape_result)
        shape_count = shape_count + 1

    return shapes_buffer


# noinspection PyPep8
# <Material diffuseColor="0.8 0.8 0.95" shininess="0.1" specularColor="0.75 0.75 0.75" ambientIntensity="0.2" emissiveColor="0.0,0.0,0.0"></Material>

# noinspection PyPep8
MODEL_TEMPLATE = """<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE X3D PUBLIC "ISO//Web3D//DTD X3D 3.3//EN" "http://www.web3d.org/specifications/x3d-3.3.dtd">
<X3D profile="Immersive" version="3.3" xsd:noNamespaceSchemaLocation="http://www.web3d.org/specifications/x3d-3.3.xsd" xmlns:xsd="http://www.w3.org/2001/XMLSchema-instance">
    <head>
        <meta content="{{file}}" name="{{title}}"/>
        <meta content="generated from {{function}}" name="{{summary}}"/>
        <meta content="{{date}}" name="created"/>
    </head>
    <Scene>{{shapes}}
    </Scene>
</X3D>
"""

# noinspection PyPep8
SHAPE_TEMPLATE = """
        <Shape DEF='{{shape_name}}'>
            <Appearance sortType="auto">
                <Material diffuseColor='{{color}}' shininess='0.5' specularColor='0.75 0.75 0.75' ambientIntensity='0.2'></Material>
            </Appearance>
            {{geometry_nodes}}
        </Shape>
"""

# noinspection PyPep8
INDEXED_TRIANGLE_SET_TEMPLATE = """<IndexedTriangleSet index='{{coord_index}}' solid='true' normalPerVertex='false'>
                <Coordinate point='{{point}}'/>
            </IndexedTriangleSet>
"""

# noinspection PyPep8
BINARY_GEOMETRY_TEMPLATE = """<BinaryGeometry DEF='BG_1' vertexCount='{{vertex_count}}' solid='true' normalPerVertex='false' primType='"TRIANGLES"' position='0.0 0.0 0.0' size='1.0 1.0 1.0' index='{{index}}' coord='{{coord}}' normal='{{normal}}'></BinaryGeometry>"""

# """<binaryGeometry solid='false' lit='false'  DEF='BG_1' vertexCount='198' primType='"TRIANGLES"' position='3977.91992188 237.727508545 -1830.51000977' size='15.1198730469 11.2089996338 18.8000488281' index='binGeo/BG_1_indexBinary.bin' coord='binGeo/BG_1_coordBinary.bin' normal='binGeo/BG_1_normalBinary.bin' texCoord='binGeo/BG_1_texCoordBinary.bin' color='binGeo/BG_1_colorBinary.bin'></binaryGeometry>"""

# # ambientIntensity='0.2' transparency='0.01'/>

# SOME example materials for future use
# diffuseColor = '0.4 0.4 0.4'
# shininess = '0.117188'
# specularColor = '0.75 0.75 0.75'
#
# diffuseColor = '1 0.949004 0.513726'
# shininess = '0.117188'
# specularColor = '0.75 0.75 0.75'
#
# diffuseColor = '1 1 0'
# shininess = '0.117188'
# specularColor = '0.75 0.75 0.75'
#
# diffuseColor = '0.6 0.6 0.6'
# shininess = '0.117188'
# specularColor = '0.75 0.75 0.75'
#
# diffuseColor = '1 0.980392 0.749004'
# shininess = '0.117188'
# specularColor = '0.75 0.75 0.75'
#
# diffuseColor = '0 0 0.752941'
# shininess = '0.273438'
# specularColor = '0.35 0.35 0.35'
#
# diffuseColor = '0.596063 0.666667 0.686259'
# shininess = '0.117188'
# specularColor = '0.75 0.75 0.75'
#
# diffuseColor = '0.298024 0.298024 0.298024'
# shininess = '0.117188'
# specularColor = '0.75 0.75 0.75'
#
# diffuseColor = '1 0.960769 0.654902'
# shininess = '0.117188'
# specularColor = '0.75 0.75 0.75'
#
# diffuseColor = '0.949004 0.949004 0.949004'
# shininess = '0.117188'
# specularColor = '0.75 0.75 0.75'
#
# diffuseColor = '0.6 1 0.792157'
# shininess = '0.117188'
# specularColor = '0.75 0.75 0.75'
#
# diffuseColor = '1 1 1'
# shininess = '0.117188'
# specularColor = '0.75 0.75 0.75'
#
# diffuseColor = '0.8 1 0.2'
# shininess = '0.117188'
# specularColor = '0.75 0.75 0.75'
#
# diffuseColor = '0.8 1 0.6'
# shininess = '0.117188'
# specularColor = '0.75 0.75 0.75'
#
# diffuseColor = '0.882338 0.960769 1'
# shininess = '0.117188'
# specularColor = '0.75 0.75 0.75'
#
# diffuseColor = '0.843137 1 0.780392'
# shininess = '0.117188'
# specularColor = '0.75 0.75 0.75'
#
# diffuseColor = '0.8 1 0.4'
# shininess = '0.117188'
# specularColor = '0.75 0.75 0.75'

# WRAPPER_TEMPLATE = """
#     <script type='text/javascript' src='http://www.x3dom.org/download/x3dom.js'> </script>
#     <link rel='stylesheet' type='text/css' href='http://www.x3dom.org/download/x3dom.css'/>
#     <h1>{{head}}</h1>
#     <h2>{{subhead}}</h2>
#     <p>{{description}}</p>
#     <x3d showLog='{{show_log}}' showStat='{{show_log}}' width='{{width}}px' height='{{height}}px'>
#         <scene>
#             <Viewpoint centerOfRotation="{{center}}"></Viewpoint>
#             <Inline nameSpaceName="{{name_space}}" mapDEFToID="true" url="{{scene_file}}" />
#         </scene>
#     </x3d>
# """

# noinspection PyPep8
WRAPPER_TEMPLATE = """
    <script type='text/javascript' src='{{x3dom_js}}'> </script>
    <link rel='stylesheet' type='text/css' href='{{x3dom_css}}'/>
    <h1>{{head}}</h1>
    <h2>{{subhead}}</h2>
    <p>{{description}}</p>
    <x3d showLog='{{show_log}}' showStat='{{show_log}}' width='{{width}}px' height='{{height}}px'>
        <scene>
            <Environment id='myEnv' SSAO='true' SSAOamount='0.3' SSAOblurDepthTreshold='1.0' SSAOradius='0.7' SSAOrandomTextureSize='4'></Environment>
            <Viewpoint centerOfRotation="{{center}}"></Viewpoint>
            <navigationinfo headlight="false"></navigationinfo>
            <background skycolor="0.8 0.8 0.8"></background>
            <directionallight def="KEY_LIGHT" ConsoleCode="0.9 0.9 1.0" direction="-0.7 -0.7 -0.3" intensity="0.5" on="true" shadowmapsize="1024" znear="-1" zfar="-1" shadowcascades="1" shadowsplitfactor="1" shadowsplitoffset="0.1"></directionallight>
            <directionallight def="FILL_LIGHT" ConsoleCode="0.9 0.7 0.4" direction="0.7  0.7 -0.3" intensity="0.4" on="true" shadowmapsize="1024" znear="-1" zfar="-1" shadowcascades="1" shadowsplitfactor="1" shadowsplitoffset="0.1"></directionallight>
            <directionallight def="BACK_LIGHT" ConsoleCode="1.0 0.9 0.0" direction="0.0  0.7  0.7" intensity="0.2" on="true" shadowmapsize="1024" znear="-1" zfar="-1" shadowcascades="1" shadowsplitfactor="1" shadowsplitoffset="0.1"></directionallight>
            <Inline nameSpaceName="{{name_space}}" mapDEFToID="true" url="{{scene_file}}" />
        </scene>
    </x3d>
"""

# noinspection PyPep8
X3D_GLTF_TEMPLATE = """
<script src="{{x3dom_js}}"></script>
<x3d id="x3d" width="800px" height="800px"> 
    <scene>
        <Environment frustumCulling='false'></Environment>
        <Viewpoint id="vp" position="0 0 5" zNear="0.01" zFar="10000"></Viewpoint>
        <navigationInfo id="head" headlight='true' type='"EXAMINE"'>  </navigationInfo> 
        <Transform DEF="model" rotation="0 1 0 0" id="gltf">
            <!-- <ExternalShape id="exshape" url="https://cdn.rawgit.com/cx20/gltf-test/master/sampleModels/Duck/glTF-Binary/Duck.glb"></ExternalShape> -->
            <!-- <ExternalShape id="exshape" url="https://cdn.rawgit.com/cx20/gltf-test/master/sampleModels/GearboxAssy/glTF-Binary/GearboxAssy.glb"></ExternalShape> -->
            <ExternalShape id="exshape" url="Box.gltf"></ExternalShape>
        </Transform>                
        <timeSensor DEF="time" cycleInterval="40" loop="true"></timeSensor>
        <OrientationInterpolator DEF="move" key="0 0.5 1" keyValue="0 1 0 0, 0 1 0 3.14, 0 1 0 6.28"></OrientationInterpolator>
        <Route fromNode="time" fromField ="fraction_changed" toNode="move" toField="set_fraction"></Route>
        <Route fromNode="move" fromField ="value_changed" toNode="model" toField="rotation"></Route>
    </scene> 
</x3d> 
"""

# noinspection PyPep8
GLTF_TEMPLATE = """
{
    "asset": {
        "generator": "affect",
        "version": "0.1"
    },
    "scene": 0,
    "scenes": [
        {
            "nodes": [
                0
            ]
        }
    ],
    "nodes": [
        {
            "children": [
                1
            ],
            "matrix": [
                1.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                -1.0,
                0.0,
                0.0,
                1.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                1.0
            ]
        },
        {
            "mesh": 0
        }
    ],
    "meshes": [
        {
            "primitives": [
                {
                    "attributes": {
                        "NORMAL": 1,
                        "POSITION": 2
                    },
                    "indices": 0,
                    "mode": 4,
                    "material": 0
                }
            ],
            "name": "Mesh"
        }
    ],
    "accessors": [
        {
            "bufferView": 0,
            "byteOffset": 0,
            "componentType": 5123,
            "count": 36,
            "max": [
                23
            ],
            "min": [
                0
            ],
            "type": "SCALAR"
        },
        {
            "bufferView": 1,
            "byteOffset": 0,
            "componentType": 5126,
            "count": 24,
            "max": [
                1.0,
                1.0,
                1.0
            ],
            "min": [
                -1.0,
                -1.0,
                -1.0
            ],
            "type": "VEC3"
        },
        {
            "bufferView": 1,
            "byteOffset": 288,
            "componentType": 5126,
            "count": 24,
            "max": [
                0.5,
                0.5,
                0.5
            ],
            "min": [
                -0.5,
                -0.5,
                -0.5
            ],
            "type": "VEC3"
        }
    ],
    "materials": [
        {
            "pbrMetallicRoughness": {
                "baseColorFactor": [
                    0.800000011920929,
                    0.0,
                    0.0,
                    1.0
                ],
                "metallicFactor": 0.0
            },
            "name": "Red"
        }
    ],
    "bufferViews": [
        {
            "buffer": 0,
            "byteOffset": 576,
            "byteLength": 72,
            "target": 34963
        },
        {
            "buffer": 0,
            "byteOffset": 0,
            "byteLength": 576,
            "byteStride": 12,
            "target": 34962
        }
    ],
    "buffers": [
        {
            "byteLength": 648,
            "uri": "Box0.bin"
        }
    ]
}
"""
