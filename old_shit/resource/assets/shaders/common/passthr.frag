#version 150

in vec4 Color;
out vec4 oColor;
void main(void)
{
    oColor = vec4( 1 ) * Color;
}