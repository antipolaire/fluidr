#version 330

layout (points) in;
layout (line_strip, max_vertices = 2) out;

in vec3 Dir[];
out vec4 Color;


vec3 computeColor(float normal_value) {
    vec3 color;
    if (normal_value < 0.0) normal_value = 0.0;
    if (normal_value > 1.0) normal_value = 1.0;
    float v1 = 1.0 / 7.0;
    float v2 = 2.0 / 7.0;
    float v3 = 3.0 / 7.0;
    float v4 = 4.0 / 7.0;
    float v5 = 5.0 / 7.0;
    float v6 = 6.0 / 7.0;
    //compute color
    if (normal_value < v1) {
        float c = normal_value / v1;
        color.x = 70. * (1. - c);
        color.y = 70. * (1. - c);
        color.z = 219. * (1. - c) + 91. * c;
    } else if (normal_value < v2) {
        float c = (normal_value - v1) / (v2 - v1);
        color.x = 0.;
        color.y = 255. * c;
        color.z = 91. * (1. - c) + 255. * c;
    } else if (normal_value < v3) {
        float c = (normal_value - v2) / (v3 - v2);
        color.x = 0. * c;
        color.y = 255. * (1. - c) + 128. * c;
        color.z = 255. * (1. - c) + 0. * c;
    } else if (normal_value < v4) {
        float c = (normal_value - v3) / (v4 - v3);
        color.x = 255. * c;
        color.y = 128. * (1. - c) + 255. * c;
        color.z = 0.;
    } else if (normal_value < v5) {
        float c = (normal_value - v4) / (v5 - v4);
        color.x = 255. * (1. - c) + 255. * c;
        color.y = 255. * (1. - c) + 96. * c;
        color.z = 0.;
    } else if (normal_value < v6) {
        float c = (normal_value - v5) / (v6 - v5);
        color.x = 255. * (1. - c) + 107. * c;
        color.y = 96. * (1. - c);
        color.z = 0.;
    } else {
        float c = (normal_value - v6) / (1. - v6);
        color.x = 107. * (1. - c) + 223. * c;
        color.y = 77. * c;
        color.z = 77. * c;
    }
    return color;
}

void main()
{
    int i;
    for (i = 0;i < gl_in.length();i++){
        vec3 dir = Dir[i];
        vec4 col = vec4(computeColor(length(dir)),1.0);
        // vec4 col = vec4(vec3(clamp(1.0, 0.0, 1.0)), 1.0);
        gl_Position = gl_in[i].gl_Position;
        Color = col;
        EmitVertex();
        gl_Position = gl_in[i].gl_Position + vec4(clamp(dir, vec3(-2.0), vec3(2.0)), 0.0) * 10.0;
        Color = col;
        EmitVertex();
        EndPrimitive();
    }
}