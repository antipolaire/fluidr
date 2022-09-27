#version 150

uniform mat4 ciModelViewProjection;

in vec4 ciPosition;
in vec3 direction;
out lowp vec3 Dir;
void main(void)
{
    gl_Position = ciModelViewProjection * ciPosition;
    Dir = direction;
}