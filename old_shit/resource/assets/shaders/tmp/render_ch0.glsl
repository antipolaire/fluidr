#version 410

uniform sampler2D iChannel0;

in highp vec2 vertTexCoord0;

uniform vec3 iResolution;

out vec4 fragColor;

void main()
{
    fragColor = texture(iChannel0,vertTexCoord0);
}
