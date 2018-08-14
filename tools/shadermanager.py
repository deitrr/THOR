"""
Manage Shader programs and associated variables

TODO: might not be the best way to manage multiple different objects.
      Need to learn more about how to structure rendering and what is optimal
"""

from PyQt5.QtGui import (QOpenGLShader,
                         QOpenGLShaderProgram,
                         QMatrix4x4,
                         QVector4D)


class ShaderManager:

    vertex_shader = """#version 400\n
                       layout(location = 0) in vec3 vp;
                       layout(location = 1) in vec3 vertex_colour;
                       layout(location = 2) in vec3 normals;

                       out vec3 normal;
                       out vec3 FragPos;
                       out vec3 colour_;
                       uniform highp mat4 view;
                       uniform highp mat4 projection;
                       uniform highp mat4 model;
                       vec3 lightdir;
                       float NdotL;
                       void main() {
                          // light dir
                          lightdir = normalize(vec3(1.0,1.0,-1.0));
                          normal = normals;
                          NdotL = 0.3+0.7*max(dot(normal, lightdir), 0.0);
                          //normal = projection * \
                          //    view * model *  vec4(normals, 1.0);
                          gl_Position = projection * \
                              view * model *  vec4(vp, 1.0);
                         FragPos = vec3(model * vec4(vp, 1.0));
                         colour_ = NdotL*vertex_colour;
                       }"""

    fragment_shader = """#version 400\n
                        in vec3 normal;
                        in vec3 FragPos;
                        in vec3 colour_;
                        out vec4 frag_colour;
                        uniform highp vec4 colour;
                        void main() {

                           frag_colour = vec4(colour_, 1.0);
                        }"""

    def __init__(self, parent):
        self.parent = parent

        # self.m_program = QOpenGLShaderProgram(parent)

        # if self.m_program.addShaderFromSourceCode(QOpenGLShader.Vertex,
        #                                           self.vertex_shader):
        #     print("initialised vertex shader")
        # else:
        #     print("vertex shader failed")

        # self.m_program.addShaderFromSourceCode(QOpenGLShader.Fragment,
        #                                        self.fragment_shader)

        # self.m_program.link()

        # self.m_viewUniform = self.m_program.uniformLocation('view')
        # self.m_projectionUniform = self.m_program.uniformLocation(
        #     'projection')
        # self.m_modelUniform = self.m_program.uniformLocation('model')

        # self.m_colourUniform = self.m_program.uniformLocation('colour')

        # # colour and position shader
        # self.m_program_pc = QOpenGLShaderProgram(parent)

        # if self.m_program_pc.addShaderFromSourceCode(QOpenGLShader.Vertex,
        #                                              self.vertex_shader_colour):
        #     print("initialised vertex color shader")
        # else:
        #     print("vertex shader color failed")

        # if self.m_program_pc.addShaderFromSourceCode(QOpenGLShader.Geometry,
        #                                              self.geometry_colour):
        #     print("initialised geometryshader")
        # else:
        #     print("geometry shader failed")
        # self.m_program_pc.addShaderFromSourceCode(QOpenGLShader.Fragment,
        #                                           self.vertex_fragment_colour)

        # self.m_program_pc.link()

        # self.m_viewUniform_pc = self.m_program_pc.uniformLocation('view')
        # self.m_projectionUniform_pc = self.m_program_pc.uniformLocation(
        #     'projection')
        # self.m_modelUniform_pc = self.m_program_pc.uniformLocation('model')
        # self.m_wirelimitUniform = self.m_program_pc.uniformLocation(
        #     'wire_limit')
        # self.m_program_field = QOpenGLShaderProgram(parent)

        # if self.m_program_field.addShaderFromSourceCode(QOpenGLShader.Vertex,
        #                                                 self.vertex_shader_field):
        #     print("initialised vertex shader")
        # else:
        #     print("vertex shader failed")

        # self.m_program_field.addShaderFromSourceCode(QOpenGLShader.Fragment,
        #                                              self.fragment_shader_field)

        # self.m_program_field.link()

        # self.m_viewUniform_field = self.m_program_field.uniformLocation('view')
        # self.m_projectionUniform_field = self.m_program_field.uniformLocation(
        #     'projection')
        # self.m_modelUniform_field = self.m_program_field.uniformLocation(
        #     'model')

        # self.m_colourUniform_field = self.m_program_field.uniformLocation(
        #     'colour')

        self.shader_dict = {}

    def add_shaders(self, shader_name, vertex=None, geometry=None, fragment=None, uniforms=[]):
        program = QOpenGLShaderProgram(self.parent)
        if vertex is not None:
            program.addShaderFromSourceCode(QOpenGLShader.Vertex, vertex)
        if geometry is not None:
            program.addShaderFromSourceCode(QOpenGLShader.Geometry, geometry)
        if fragment is not None:
            program.addShaderFromSourceCode(QOpenGLShader.Fragment, fragment)

        program.link()

        uniforms_dict = {}
        for u in uniforms:
            uniforms_dict[u] = program.uniformLocation(u)

        self.shader_dict[shader_name] = {'program': program,
                                         'uniforms': uniforms_dict}

    def load(self, name):
        self.shader_dict[name]['program'].bind()

    def release(self, name):
        self.shader_dict[name]['program'].release()

    def set_uniform(self, shader_name, uniform_name, uniform_value):
        program = self.shader_dict[shader_name]['program']
        uniform = self.shader_dict[shader_name]['uniforms'][uniform_name]
        program.setUniformValue(uniform, uniform_value)

    # def start_paint_field(self):
    #     self.m_program_field.bind()

    # def end_paint_field(self):
    #     self.m_program_field.release()

    # def start_paint(self):
    #     self.m_program.bind()

    # def end_paint(self):
    #     self.m_program.release()

    # def start_paint_pc(self):
    #     self.m_program_pc.bind()

    # def end_paint_pc(self):
    #     self.m_program_pc.release()

    # def set_colour(self, colour):
    #     self.m_program.setUniformValue(self.m_colourUniform, colour)

    # def set_colour_field(self, colour):
    #     self.m_program_field.setUniformValue(
    #         self.m_colourUniform_field, colour)

    # def set_wirelimit(self, wirelimit):
    #     self.m_program_pc.setUniformValue(self.m_wirelimitUniform, wirelimit)

    # def set_view(self, view):
    #     self.m_program.setUniformValue(
    #         self.m_viewUniform, view)

    # def set_view_pc(self, view):
    #     self.m_program_pc.setUniformValue(
    #         self.m_viewUniform_pc, view)

    # def set_view_field(self, view):
    #     self.m_program_field.setUniformValue(
    #         self.m_viewUniform_field, view)

    # def set_projection(self, projection):
    #     self.m_program.setUniformValue(
    #         self.m_projectionUniform, projection)

    # def set_projection_field(self, projection):
    #     self.m_program_field.setUniformValue(
    #         self.m_projectionUniform_field, projection)

    # def set_projection_pc(self, projection):
    #     self.m_program_pc.setUniformValue(
    #         self.m_projectionUniform_pc, projection)

    # def set_model(self, model):
    #     self.m_program.setUniformValue(
    #         self.m_modelUniform, model)

    # def set_model_field(self, model):
    #     self.m_program_field.setUniformValue(
    #         self.m_modelUniform_field, model)

    # def set_model_pc(self, model):
    #     self.m_program_pc.setUniformValue(
    #         self.m_modelUniform_pc, model)
