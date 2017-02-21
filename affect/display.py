from jinja2 import Template
from IPython.core.display import HTML

wrapper_template = """
    <script type='text/javascript' src='http://www.x3dom.org/download/x3dom.js'> </script>
    <link rel='stylesheet' type='text/css' href='http://www.x3dom.org/download/x3dom.css'/>
    <h1>{{head}}</h1>
    <h2>{{subhead}}</h2>
    <p>{{description}}</p>
    <x3d width='{{width}}px' height='{{height}}px'>
        <scene>
            <navigationinfo headlight="false"></navigationinfo>
            <background skycolor="0.8 0.8 0.8"></background>
            <directionallight def="KEY_LIGHT" ConsoleCode="0.9 0.9 1.0" direction="-0.7 -0.7 -0.3" intensity="1.0" on="true" shadowmapsize="1024" znear="-1" zfar="-1" shadowcascades="1" shadowsplitfactor="1" shadowsplitoffset="0.1"></directionallight>
            <directionallight def="FILL_LIGHT" ConsoleCode="0.5 0.7 0.8" direction="0.7  0.7 -0.3" intensity="0.9" on="true" shadowmapsize="1024" znear="-1" zfar="-1" shadowcascades="1" shadowsplitfactor="1" shadowsplitoffset="0.1"></directionallight>
            <directionallight def="BACK_LIGHT" ConsoleCode="0.7 0.9 1.0" direction="0.0  0.7  0.7" intensity="0.4" on="true" shadowmapsize="1024" znear="-1" zfar="-1" shadowcascades="1" shadowsplitfactor="1" shadowsplitoffset="0.1"></directionallight>
            <Inline nameSpaceName="Model" mapDEFToID="true" url="{{scene_file}}" />
        </scene>
    </x3d>
"""

def notebook_scene(head, subhead, description, width, height, scene_file):
    """
    Return valid HTML code that will display interactive 3D model using the X3DOM framework.

    An example of displaying a model inside an ipython notebook would be as follows.

    :Example:

            >>> from affect import display.notebook_scene
            >>> notebook_scene('Model', 'A solid mechanics model', 'analysis series 3c', 500, 500, 'myModel.x3d')

    :param my_head: a string with a title for the model
    :param my_subhead: a string with a subtitle for the model
    :param my_description: longer text description for this model view
    :param width: width of the display window in pixels
    :param height: height of the display window in pixels
    :param scene_file: an existing x3d scene using the schema at "http://www.web3d.org/specifications/x3d-3.1.xsd'
    :return: string -- HTML code useful for embedding inside a <body> tag or as part of a ipython notebook.
    """
    t = Template(wrapper_template)
    source = t.render(head=head, subhead=subhead, description=description, width=width, height=height,
                      scene_file=scene_file)
    return HTML(source)

model_template = """<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE X3D PUBLIC "ISO//Web3D//DTD X3D 3.1//EN" "http://www.web3d.org/specifications/x3d-3.1.dtd">
<X3D profile="Immersive" version="3.1" xsd:noNamespaceSchemaLocation="http://www.web3d.org/specifications/x3d-3.1.xsd" xmlns:xsd="http://www.w3.org/2001/XMLSchema-instance">
    <head>
        <meta content="{{file}}" name="{{title}}"/>
        <meta content="generated from {{function}}" name="{{summary}}"/>
        <meta content="{{date}}" name="created"/>
    </head>
    <Scene>
        <Shape>
            <Appearance sortType="auto">
                <Material diffuseColor="0.8 0.8 0.95" shininess="0.1" specularColor="0.75 0.75 0.75" ambientIntensity="0.2" emissiveColor="0.0,0.0,0.0"></Material>
            </Appearance>
            <IndexedFaceSet coordIndex='0 1 5 4 -1 1 2 6 5 -1 2 3 7 6 -1 3 0 4 7 -1 0 3 2 1 -1 4 5 6 7 -1'>
                <Coordinate point='-1.0 -1.0 1.0 1.0 -1.0 1.0 1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 -1.0 -1.0 1.0 -1.0'/>
            </IndexedFaceSet>
        </Shape>
    </Scene>
</X3D>
"""

def write_scene(file, title, summary, faces, coordinates):
    import time
    func_name = 'affect.display.write_scene'
    now = time.strftime("%H:%M:%S")
    date = time.strftime("%Y-%m-%d")
    t = Template(model_template)
    result = t.render(file=file, title=title, function=func_name, summary=summary, date=date, time=now)
    with open(file, "w") as text_file:
        text_file.write(result)
