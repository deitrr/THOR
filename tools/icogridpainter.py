
from PyQt5.QtGui import (QOpenGLShader,
                         QOpenGLShaderProgram,
                         QMatrix4x4,
                         QVector4D,
                         QPainterPath,
                         QFont,
                         QSurfaceFormat)

import OpenGL.GL as gl
import numpy as np
import ctypes
import math
from math import cos, sin, pi, pow, sqrt

import matplotlib.pyplot as plt

from utilities import HSV_to_RGB, spherical

from basepainter import BasePainter

vertex_shader_icos = """#version 400\n
                       layout(location = 0) in vec3 vp;
                       layout(location = 1) in vec3 vertex_colour;
                       layout(location = 2) in float vertex_value;

                       out VS_OUT {
                           out vec3 colour;
                           out float val;
                       } vs_out;


                       uniform highp mat4 view;
                       uniform highp mat4 projection;
                       uniform highp mat4 model;
                       void main() {
                         vs_out.colour = vertex_colour;
                         vs_out.val = vertex_value;
                          //colour = vec3(1.0,1.0,0.0);
                          gl_Position = projection * \
                              view * model *  vec4(vp, 1.0);
                       }"""

geometry_shader_icos = """
#version 400\n


layout (triangles) in;
layout (triangle_strip, max_vertices = 3) out;

in VS_OUT {
   vec3 colour;
   float val;
} gs_in[];

out GS_OUT {
   vec3 fcolor;
   vec3 dist;
   float val;
} gs_out;

void main(void)
{
   float WIN_SCALE = 3200.0;
   // taken from 'Single-Pass Wireframe Rendering'
   vec2 p0 = WIN_SCALE * gl_in[0].gl_Position.xy/gl_in[0].gl_Position.w;
   vec2 p1 = WIN_SCALE * gl_in[1].gl_Position.xy/gl_in[1].gl_Position.w;
   vec2 p2 = WIN_SCALE * gl_in[2].gl_Position.xy/gl_in[2].gl_Position.w;
   
   vec2 v0 = p2-p1;
   vec2 v1 = p2-p0;
   vec2 v2 = p1-p0;

   float area = abs(v1.x*v2.y - v1.y * v2.x);

   gs_out.fcolor = gs_in[0].colour;

   //gs_out.dist = vec3(area/length(v0),0,0);
   gs_out.dist = vec3(1.0,0,0);
   gs_out.val = gs_in[0].val;
   gl_Position = 1.5*gl_in[0].gl_Position;
   
   EmitVertex();
   gs_out.fcolor = gs_in[1].colour;
//   gs_out.dist = vec3(0,area/length(v1),0);
   gs_out.dist = vec3(0,1.0,0);
   gs_out.val = gs_in[1].val;
   gl_Position = gl_in[1].gl_Position;
   EmitVertex();
   gs_out.fcolor = gs_in[2].colour;
//   gs_out.dist = vec3(0,0,area/length(v2));
   gs_out.dist = vec3(0,0,1.0);
   gs_out.val = gs_in[2].val;
   gl_Position = gl_in[2].gl_Position;
   EmitVertex();
   EndPrimitive();
}"""

fragment_shader_icos = """
#version 400\n
in GS_OUT {
   vec3 fcolor;
   vec3 dist;
   float val;
} fs_in;

vec3 colormap(float H)
{
   vec3 RGB;
   RGB.x = H;
   RGB.y = H;
   RGB.z = H;

   return RGB;
} 

vec3 H_to_RGB(float H)
{
   H *= 360.0;
   float S = 1.0;
   float V = 1.0; 

   float C = V*S;
   float H_p = H/60.0;
   float X = C*(1.0 - abs(mod(H_p,2.0) - 1.0));

   vec3 R1G1B1;
   R1G1B1.x = 0.0;
   R1G1B1.y = 0.0;
   R1G1B1.z = 0.0;
        
    if (0 <= H_p && H_p <= 1.0)
    {
        R1G1B1.x = C;
        R1G1B1.y = X;
        R1G1B1.z = 0.0;
    }
    else if (1.0 <= H_p && H_p <= 2.0)
    {
        R1G1B1.x = X;
        R1G1B1.y = C;
        R1G1B1.z = 0.0;
    }
    else if (2.0 <= H_p && H_p <= 3.0)
    {
        R1G1B1.x = 0.0;
        R1G1B1.y = C;
        R1G1B1.z = X;
    }
    else if (3. <= H_p && H_p <= 4.0)
    {
        R1G1B1.x = 0.0;
        R1G1B1.y = X;
        R1G1B1.z = C;
    }
    else if (4.0 <= H_p && H_p <= 5.0)
    {
        R1G1B1.x = X;
        R1G1B1.y = 0.0;
        R1G1B1.z = C;
    }
    else if (5.0 <= H_p && H_p < 6.0)
    {
        R1G1B1.x = C;
        R1G1B1.y = 0.0;
        R1G1B1.z = X;
    }
    
    float m = V-C;
    

    vec3 RGB;
    RGB.x = R1G1B1.x + m;
    RGB.y = R1G1B1.y + m;
    RGB.z = R1G1B1.z + m;
    
    return RGB;
}

out vec4 frag_colour;
uniform float wire_limit;
void main() {
   float nearD = min(min(fs_in.dist[0],fs_in.dist[1]),fs_in.dist[2]);
   float edgeIntensity = exp2(-1.0*nearD*nearD);
   if (nearD < wire_limit)
      frag_colour = vec4(0.1, 0.1, 0.1, 1.0 );
   else
        //frag_colour = vec4(H_to_RGB(fs_in.val), 1.0);
//      frag_colour = vec4(colormap(fs_in.val), 1.0);
      frag_colour = vec4(fs_in.fcolor, 1.0);
      //frag_colour = vec4(1.0,0.0,0.0, 1.0);
      //                                   frag_colour = (edgeIntensity * vec4( 0.1, 0.1, 0.1, 1.0 )) +
      //                                                 (1.0-edgeIntensity)*vec4(colour, 1.0);
}"""


class IcoGridPainter(BasePainter):
    def __init__(self):
        BasePainter.__init__(self)
        self.level = 0
        self.wireframe = False
        self.draw_idx = 0
        self.altitude = 0
        self.show_data = True

    def set_grid_data(self, dataset):
        self.dataset = dataset

    def set_shader_manager(self, shader_manager):
        BasePainter.set_shader_manager(self, shader_manager)

        self.shader_manager.add_shaders("icos_grid",
                                        vertex=vertex_shader_icos,
                                        geometry=geometry_shader_icos,
                                        fragment=fragment_shader_icos,
                                        uniforms=["model", "view", "projection", "wire_limit"])

    def initializeGL(self):

        # create

        # num_vertices = 12

        # vertices = np.zeros((num_vertices, 3),
        #                     dtype=np.float32)
        lonlat = self.dataset.get_grid()
        num_points, num_levels = self.dataset.get_dim()

        vertices = np.zeros((num_points, 3),
                            dtype=np.float32)
        print("number of points: ", num_points)
        # r = 1.2
        # w = 2.0*math.acos(1.0/(2.0*sin(pi/5.0)))

        # angles = [
        #     [0.0, pi/2.0],
        #     [0.0, -pi/2.0],
        #     [-pi/5.0, pi/2.0-w],
        #     [pi/5.0, pi/2.0-w],
        #     [3.0*pi/5.0, pi/2.0-w],
        #     [pi, pi/2.0-w],
        #     [-3.0*pi/5.0, pi/2.0-w],
        #     [0.0, -(pi/2.0-w)],
        #     [2.0*pi/5.0, -(pi/2.0-w)],
        #     [4.0*pi/5.0, -(pi/2.0-w)],
        #     [-4.0*pi/5.0, -(pi/2.0-w)],
        #     [-2.0*pi/5.0, -(pi/2.0-w)]
        # ]

        # indexing function through rhombis
        # level
        g = int(pow((num_points - 2)/10, 1/4)) - 2
        num_rhombi = 10

        num_points_side_region = int(pow(2, 4))
        nl_reg = num_points_side_region
        nl2 = int(pow(num_points_side_region, 2))
        # kxl = int(sqrt(num_subrhombi))
        kxl = int(pow((num_points - 2)/10, 1/2))//num_points_side_region

        def idx(fc, kx, ky, i, j):
            return nl2*(fc*nfaces + ky*kxl + kx) + j*nl_reg + i

        # nfaces
        num_subrhombi = kxl*kxl
        nfaces = num_subrhombi

        print("subdivision level: ", g)
        print("number of subregions square:", kxl)
        print("num sub rhombis", num_subrhombi)

        # build triangles
        # triangles in rhombi and subrhombi
        num_triangles = int(pow(nl_reg - 1, 2))*2*num_subrhombi*num_rhombi
        # halos
        # each subrhombus takes care of two of its borders,
        # nl_reg-1 squares, * 2 for triangles, * 2 for two borders

        # neighbours of major rhombi face indexes
        # use clockwise addressing, top left, top right, bottom right, bottom left
        rhombi_neighbours = [
            [4, 1, 5, 9],  # idx: 0
            [0, 2, 6, 5],  # idx: 1
            [1, 3, 7, 6],  # idx: 2
            [2, 4, 8, 7],  # idx: 3
            [3, 0, 9, 8],  # idx: 4

            [0, 1, 6, 9],  # idx: 5
            [1, 2, 7, 5],  # idx: 6
            [2, 3, 8, 6],  # idx: 7
            [3, 4, 9, 7],  # idx: 8
            [4, 0, 5, 8],  # idx: 9
        ]
        num_triangles += 2*(nl_reg-1)*2*num_subrhombi*num_rhombi
        # corner between two halos
        num_triangles += 2*num_subrhombi*num_rhombi
        # poles
        num_triangles += 2*5
        print("size of triangles buffer", num_triangles)

        triangles = np.zeros((num_triangles, 3), dtype=np.uint32)

        draw_rhomb = [0, 1, 2, 3, 4, 5, 6, 7, 8,  9]
        # upper
        # draw_rhomb = [0, 1, 2, 3, 4]
        # lower
        #draw_rhomb = [5, 6, 7, 8,  9]
        # North South link
        # draw_rhomb = [0, 1, 5, 6, 9]
        faces = True

        # halos = False
        halos = True
        poles = True
        corners = True
        triangle_idx = 0

        # index for corner point
        cn = nl_reg - 1

        for fc in range(num_rhombi):
            if fc not in draw_rhomb:
                continue
            # sub rhombis
            for kx in range(kxl):
                for ky in range(kxl):
                    if faces:
                        # inside one sub-rombi
                        # loop on triangles and add them, two per "square"
                        for i in range(nl_reg-1):
                            for j in range(nl_reg-1):
                                i1 = idx(fc, kx, ky, i, j)
                                i2 = idx(fc, kx, ky, i, j+1)
                                i3 = idx(fc, kx, ky, i+1, j)
                                i4 = idx(fc, kx, ky, i+1, j+1)
                                triangles[triangle_idx][0] = i1
                                triangles[triangle_idx][1] = i2
                                triangles[triangle_idx][2] = i4

                                triangle_idx += 1
                                triangles[triangle_idx][0] = i1
                                triangles[triangle_idx][1] = i4
                                triangles[triangle_idx][2] = i3
                                triangle_idx += 1

                    if corners:
                        # corner indexes
                        fc_t_l = rhombi_neighbours[fc][0]
                        fc_b_l = rhombi_neighbours[fc][3]

                        if fc < 5:
                            # upper rhombi
                            if kx == 0 and ky == 0:
                                i_c = idx(fc, 0, 0, 0, 0)
                                i_c_t = idx(fc_b_l, 0, kxl-1, 0, cn)
                                i_c_b = idx(fc_t_l,
                                            kxl - 1, kxl - 1, cn, cn)

                                triangles[triangle_idx][0] = i_c
                                triangles[triangle_idx][1] = i_c_t
                                triangles[triangle_idx][2] = i_c_b
                                triangle_idx += 1

                            elif kx == 0:
                                # on ky = 0 side
                                i_c = idx(fc, kx, ky, 0, 0)
                                i_c_t = idx(fc, kx, ky-1, 0, cn)
                                i_c_b = idx(fc_t_l,
                                            kxl - 1 - ky + 1, kxl - 1, 0, cn)

                                triangles[triangle_idx][0] = i_c
                                triangles[triangle_idx][1] = i_c_t
                                triangles[triangle_idx][2] = i_c_b
                                triangle_idx += 1

                                i_c = idx(fc, kx, ky, 0, 0)
                                i_c_t = idx(fc_t_l, kxl - 1 - ky + 1, kxl-1,
                                            0, cn)
                                i_c_b = idx(fc_t_l, kxl - 1 - ky, kxl-1,
                                            cn, cn)

                                triangles[triangle_idx][0] = i_c
                                triangles[triangle_idx][1] = i_c_t
                                triangles[triangle_idx][2] = i_c_b
                                triangle_idx += 1
                            elif ky == 0:
                                # on ky = 0 side
                                i_c = idx(fc, kx, ky, 0, 0)
                                i_c_t = idx(fc, kx-1, ky, cn, 0)
                                i_c_b = idx(fc_b_l, kx-1, kxl-1,
                                            cn, cn)

                                triangles[triangle_idx][0] = i_c
                                triangles[triangle_idx][1] = i_c_t
                                triangles[triangle_idx][2] = i_c_b
                                triangle_idx += 1

                                i_c = idx(fc, kx, ky, 0, 0)

                                i_c_t = idx(fc_b_l, kx, kxl-1,
                                            0, cn)
                                i_c_b = idx(fc_b_l, kx-1, kxl-1,
                                            cn, cn)

                                triangles[triangle_idx][0] = i_c
                                triangles[triangle_idx][1] = i_c_t
                                triangles[triangle_idx][2] = i_c_b
                                triangle_idx += 1
                            else:  # kx >0 and ky > 0
                                # inside a face
                                i_c = idx(fc, kx, ky, 0, 0)
                                i_c_t = idx(fc, kx-1, ky-1, cn, cn)
                                i_c_b = idx(fc, kx, ky-1, 0, cn)
                                triangles[triangle_idx][0] = i_c
                                triangles[triangle_idx][1] = i_c_t
                                triangles[triangle_idx][2] = i_c_b
                                triangle_idx += 1
                                i_c = idx(fc, kx, ky, 0, 0)
                                i_c_t = idx(fc, kx-1, ky, cn, 0)
                                i_c_b = idx(fc, kx-1, ky-1, cn, cn)
                                triangles[triangle_idx][0] = i_c
                                triangles[triangle_idx][1] = i_c_t
                                triangles[triangle_idx][2] = i_c_b
                                triangle_idx += 1
                        else:
                            # lower rhombi
                            if kx == 0 and ky == 0:
                                i_c = idx(fc, 0, 0, 0, 0)
                                i_c_t = idx(fc_t_l, kxl-1, 0, cn, 0)
                                i_c_b = idx(fc_b_l,
                                            kxl - 1, kxl - 1, cn, cn)

                                triangles[triangle_idx][0] = i_c
                                triangles[triangle_idx][1] = i_c_t
                                triangles[triangle_idx][2] = i_c_b
                                triangle_idx += 1

                            elif kx == 0:
                                # on kx = 0 side
                                i_c = idx(fc, kx, ky, 0, 0)
                                i_c_b = idx(fc_t_l, kxl-1, ky-1,
                                            cn, cn)

                                i_c_t = idx(fc, kx, ky-1, 0, cn)
                                # i_c_b = idx(fc_t_l, kx-1, kxl-1,
                                #             cn, cn)

                                triangles[triangle_idx][0] = i_c
                                triangles[triangle_idx][1] = i_c_t
                                triangles[triangle_idx][2] = i_c_b
                                triangle_idx += 1

                                i_c = idx(fc, kx, ky, 0, 0)

                                i_c_t = idx(fc_t_l, kxl-1, ky,
                                            cn, 0)
                                i_c_b = idx(fc_t_l, kxl-1, ky-1,
                                            cn, cn)

                                triangles[triangle_idx][0] = i_c
                                triangles[triangle_idx][1] = i_c_t
                                triangles[triangle_idx][2] = i_c_b
                                triangle_idx += 1
                            elif ky == 0:
                                # on ky = 0 side
                                i_c = idx(fc, kx, ky, 0, 0)
                                i_c_t = idx(fc, kx-1, ky, cn, 0)
                                i_c_b = idx(fc_b_l, kxl-1,
                                            kxl - 1 - kx + 1, cn, 0)

                                triangles[triangle_idx][0] = i_c
                                triangles[triangle_idx][1] = i_c_t
                                triangles[triangle_idx][2] = i_c_b
                                triangle_idx += 1

                                i_c = idx(fc, kx, ky, 0, 0)
                                i_c_t = idx(fc_b_l, kxl-1,
                                            kxl - 1 - kx + 1, cn, 0)
                                i_c_b = idx(fc_b_l, kxl-1,
                                            kxl - 1 - kx, cn, cn)

                                triangles[triangle_idx][0] = i_c
                                triangles[triangle_idx][1] = i_c_t
                                triangles[triangle_idx][2] = i_c_b
                                triangle_idx += 1

                            else:  # kx >0 and ky > 0
                                # inside a face
                                i_c = idx(fc, kx, ky, 0, 0)
                                i_c_t = idx(fc, kx-1, ky-1, cn, cn)
                                i_c_b = idx(fc, kx, ky-1, 0, cn)
                                triangles[triangle_idx][0] = i_c
                                triangles[triangle_idx][1] = i_c_t
                                triangles[triangle_idx][2] = i_c_b
                                triangle_idx += 1
                                i_c = idx(fc, kx, ky, 0, 0)
                                i_c_t = idx(fc, kx-1, ky, cn, 0)
                                i_c_b = idx(fc, kx-1, ky-1, cn, cn)
                                triangles[triangle_idx][0] = i_c
                                triangles[triangle_idx][1] = i_c_t
                                triangles[triangle_idx][2] = i_c_b
                                triangle_idx += 1
                                fc_t_l = rhombi_neighbours[fc][0]
                                fc_b_l = rhombi_neighbours[fc][3]

                                i_c_b = idx(fc_t_l, kxl-1, kxl-1, nl_reg-1, 0)
                                i_c_t = idx(fc_b_l, kxl-1, kxl -
                                            1, nl_reg-1, nl_reg-1)

                    if poles:
                        # poles indexes
                        if fc < 5:
                            if kx == 0 and ky == kxl - 1:
                                i_p = vertices.shape[0] - 2
                                i_c1 = idx(fc, 0, kxl-1, 0, nl_reg-1)

                                fc_t_l = rhombi_neighbours[fc][0]
                                i_c2 = idx(fc_t_l, 0, kxl-1, 0, nl_reg-1)
                                triangles[triangle_idx][0] = i_p
                                triangles[triangle_idx][1] = i_c1
                                triangles[triangle_idx][2] = i_c2
                                triangle_idx += 1
                        else:
                            if kx == kxl - 1 and ky == 0:
                                i_p = vertices.shape[0] - 1
                                i_c1 = idx(fc, kxl-1, 0, nl_reg-1, 0)

                                fc_b_l = rhombi_neighbours[fc][3]
                                i_c2 = idx(fc_b_l, kxl-1, 0, nl_reg-1, 0)
                                triangles[triangle_idx][0] = i_p
                                triangles[triangle_idx][1] = i_c1
                                triangles[triangle_idx][2] = i_c2
                                triangle_idx += 1

                    if halos:
                        # sub rhombi halo
                        # top left halo
                        fc_top = fc
                        kx_top = kx - 1
                        if kx == 0 and fc < 5:
                            fc_top = rhombi_neighbours[fc][0]
                            kx_top = kxl - 1
                            print(kx, ky, kx_top)

                            for i in range(nl_reg-1):
                                i1 = idx(fc, kx, ky, 0, i)
                                i2 = idx(fc, kx, ky, 0, i+1)
                                i3 = idx(fc_top, kxl-1-ky, kxl-1,
                                         nl_reg - 1 - (i), nl_reg-1)
                                i4 = idx(fc_top, kxl-1-ky, kxl-1,
                                         nl_reg - 1 - (i+1), nl_reg-1)
                                triangles[triangle_idx][0] = i1
                                triangles[triangle_idx][1] = i3
                                triangles[triangle_idx][2] = i2

                                triangle_idx += 1
                                triangles[triangle_idx][0] = i2
                                triangles[triangle_idx][1] = i3
                                triangles[triangle_idx][2] = i4
                                triangle_idx += 1
                        elif kx == 0 and fc > 4:
                            fc_top = rhombi_neighbours[fc][0]
                            for i in range(nl_reg-1):
                                i1 = idx(fc, kx, ky, 0, i)
                                i2 = idx(fc, kx, ky, 0, i+1)
                                i3 = idx(fc_top, kxl-1, ky, nl_reg - 1, i)

                                i4 = idx(fc_top, kxl-1, ky, nl_reg - 1, i+1)

                                triangles[triangle_idx][0] = i1
                                triangles[triangle_idx][1] = i3
                                triangles[triangle_idx][2] = i2

                                triangle_idx += 1
                                triangles[triangle_idx][0] = i2
                                triangles[triangle_idx][1] = i3
                                triangles[triangle_idx][2] = i4
                                triangle_idx += 1
                        else:
                            for i in range(nl_reg-1):
                                i1 = idx(fc, kx, ky, 0, i)
                                i2 = idx(fc, kx, ky, 0, i+1)
                                i3 = idx(fc, kx-1, ky, nl_reg-1, i)
                                i4 = idx(fc, kx-1, ky, nl_reg-1, i+1)

                                triangles[triangle_idx][0] = i1
                                triangles[triangle_idx][1] = i3
                                triangles[triangle_idx][2] = i2

                                triangle_idx += 1
                                triangles[triangle_idx][0] = i2
                                triangles[triangle_idx][1] = i3
                                triangles[triangle_idx][2] = i4
                                triangle_idx += 1
                        if ky == 0 and fc > 4:
                            fc_bot = rhombi_neighbours[fc][3]
                            for i in range(nl_reg-1):
                                i1 = idx(fc, kx, ky,  i, 0)
                                i2 = idx(fc, kx, ky,  i+1, 0)
                                i3 = idx(fc_bot, kxl-1, kxl-1-kx,
                                         nl_reg-1, nl_reg-1 - i)
                                i4 = idx(fc_bot, kxl-1, kxl-1-kx,
                                         nl_reg-1, nl_reg-1-(i+1))
                                triangles[triangle_idx][0] = i1
                                triangles[triangle_idx][1] = i3
                                triangles[triangle_idx][2] = i2

                                triangle_idx += 1
                                triangles[triangle_idx][0] = i2
                                triangles[triangle_idx][1] = i3
                                triangles[triangle_idx][2] = i4
                                triangle_idx += 1
                        elif ky == 0 and fc < 5:
                            fc_bot = rhombi_neighbours[fc][3]
                            ky_bot = kxl - 1
                            print(kx, ky, ky_bot)
                            for i in range(nl_reg-1):
                                i1 = idx(fc, kx, ky,  i, 0)
                                i2 = idx(fc, kx, ky,  i+1, 0)
                                i3 = idx(fc_bot, kx, ky_bot, i, nl_reg-1)
                                i4 = idx(fc_bot, kx, ky_bot, i+1, nl_reg-1)
                                triangles[triangle_idx][0] = i1
                                triangles[triangle_idx][1] = i3
                                triangles[triangle_idx][2] = i2

                                triangle_idx += 1
                                triangles[triangle_idx][0] = i2
                                triangles[triangle_idx][1] = i3
                                triangles[triangle_idx][2] = i4
                                triangle_idx += 1
                        else:
                            for i in range(nl_reg-1):
                                i1 = idx(fc, kx, ky,  i, 0)
                                i2 = idx(fc, kx, ky,  i+1, 0)
                                i3 = idx(fc, kx, ky-1, i, nl_reg-1)
                                i4 = idx(fc, kx, ky-1, i+1, nl_reg-1)
                                triangles[triangle_idx][0] = i1
                                triangles[triangle_idx][1] = i3
                                triangles[triangle_idx][2] = i2

                                triangle_idx += 1
                                triangles[triangle_idx][0] = i2
                                triangles[triangle_idx][1] = i3
                                triangles[triangle_idx][2] = i4
                                triangle_idx += 1

        print("number of created triangles:", triangle_idx)

        r = 1.0
        for i in range(vertices.shape[0]):
            # spherical -> theta vertical mvt, declination-> latitude
            # -> phi -> horizontal mvt, azimuth -> longitude
            vertices[i, :] = spherical(r,
                                       lonlat[i, 1],
                                       lonlat[i, 0])

        # build VAOs for all input data
        self.vao_list = []

        # self.grid_elements_count = elements.size
        # self.grid_elements_count = lines.size
        self.grid_elements_count = triangles.size

        # self.vao_grid = gl.glGenVertexArrays(1)
        # gl.glBindVertexArray(self.vao_grid)
        # gl.glEnableVertexAttribArray(0)

        self.grid_vertex_count = vertices.size

        all_vertices = np.zeros((num_levels,  num_points, 3),
                                dtype=np.float32)
        self.all_colors = np.zeros((num_levels,  num_points, 3),
                                   dtype=np.float32)
        self.grid_scalar_data = np.zeros((num_levels,  num_points),
                                         dtype=np.float32)

        all_triangles = np.zeros((num_levels,  triangles.shape[0], 3),
                                 dtype=np.uint32)

        self.grid_color_data = np.zeros((num_levels,  num_points, 3),
                                        dtype=np.float32)
        print("num_levels:", num_levels)
        relative_radius = self.dataset.get_relative_radii()

        for fc in range(num_rhombi):
            for kx in range(kxl):
                for ky in range(kxl):
                    for i in range(nl_reg):
                        for j in range(nl_reg):
                            c = idx(fc, kx, ky,  i, j)
                            self.grid_color_data[:, c, :] = HSV_to_RGB(
                                360.0*c/(num_points - 1), 1.0, 1.0)
                            self.grid_scalar_data[:, c] = c/(num_points - 1)
#        for i in range(num_points):
#            self.grid_color_data[:, i, :] = HSV_to_RGB(
#                360.0*i/(num_points - 1), 1.0, 1.0)

        for level in range(num_levels-1):
            r = relative_radius[level]

            all_vertices[level, :, :] = r*vertices

            all_triangles[level, :, :] = triangles + \
                level*vertices.shape[0]

        self.vao_list = self.create_sphere_vao(all_vertices,
                                               all_triangles,
                                               self.all_colors,
                                               self.grid_scalar_data)

        self.num_levels = num_levels
        self.num_points = num_points
        self.last_loaded = -1

    def update_colors(self):
        if not self.show_data:

            self.update_colors_vbo(self.grid_color_data)
            self.update_values_vbo(self.grid_scalar_data)
        else:
            if self.draw_idx == self.last_loaded:
                return
            data, data_scalar = self.dataset.get_color_data(self.draw_idx)

            data = np.array(data, copy=True)
            data_scalar = np.array(data_scalar, copy=True)
            self.update_colors_vbo(data)
            self.update_values_vbo(data_scalar)
            self.last_loaded = self.draw_idx

    def create_sphere_vao(self, vertices, triangles, colors, values):
        vao = gl.glGenVertexArrays(1)

        gl.glBindVertexArray(vao)
        gl.glEnableVertexAttribArray(0)

        self.create_vbo(vertices)

        self.color_vbo = self.create_colors_vbo(colors, dynamic=True)
        self.values_vbo = self.create_values_vbo(values, dynamic=True)
        self.create_elements_vbo(triangles)

        return vao

    def paint_grid(self):
        # display grid
        if self.wireframe:
            self.shader_manager.set_uniform("icos_grid", "wire_limit", 0.00001)
        else:
            self.shader_manager.set_uniform("icos_grid", "wire_limit", -10.0)
            #            gl.glPolygonMode(gl.GL_FRONT_AND_BACK, gl.GL_LINE)
        #        print("paint grid",
        #              self.vao_list[self.draw_idx], self.grid_vertex_count)
        self.update_colors()

        gl.glBindVertexArray(self.vao_list)
        idx = 4*self.grid_elements_count * self.altitude

        # gl.glDrawElements(gl.GL_TRIANGLES,
        #                  int(self.grid_elements_count),
        #                  gl.GL_UNSIGNED_INT,
        #                  ctypes.c_void_p(idx))
        gl.glDrawElements(gl.GL_TRIANGLES,
                          int(self.grid_elements_count),
                          gl.GL_UNSIGNED_INT,
                          ctypes.c_void_p(idx))
        # if self.draw_idx == 0:
        #     gl.glBindVertexArray(self.vao_grid)
        #     gl.glDrawElements(gl.GL_TRIANGLES,
        #                       self.grid_elements_count,
        #                       gl.GL_UNSIGNED_INT,
        #                       None)
        # else:
        #     gl.glBindVertexArray(self.vao_list)

        #     idx = 4*self.grid_elements_count *
        #         (self.num_levels*(self.draw_idx-1) + self.altitude)

        #     gl.glDrawElements(gl.GL_TRIANGLES,
        #                       int(self.grid_elements_count),
        #                       gl.GL_UNSIGNED_INT,
        #                       ctypes.c_void_p(idx))
#        if self.wireframe:
#            gl.glPolygonMode(gl.GL_FRONT_AND_BACK, gl.GL_FILL)

    def update_colors_vbo(self, colors):
        # print("update colors")
        self.update_vbo(self.vao_list,
                        self.color_vbo,
                        colors)

    def update_values_vbo(self, values):
        # print("update colors")
        self.update_vbo(self.vao_list,
                        self.values_vbo,
                        values)

    def create_values_vbo(self, values, dynamic=False):
        vbo = gl.glGenBuffers(1)

        gl.glBindBuffer(gl.GL_ARRAY_BUFFER, vbo)

        if dynamic:
            gl.glBufferData(gl.GL_ARRAY_BUFFER,
                            values.nbytes,
                            values,
                            gl.GL_DYNAMIC_DRAW)
        else:
            gl.glBufferData(gl.GL_ARRAY_BUFFER,
                            values.nbytes,
                            values,
                            gl.GL_STATIC_DRAW)

        gl.glVertexAttribPointer(2, 1,
                                 gl.GL_FLOAT, False,
                                 0, None)
        gl.glEnableVertexAttribArray(2)

        return vbo


if __name__ == '__main__':
    ico = ico(4)
