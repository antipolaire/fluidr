#version 410

uniform sampler2D iChannel0;

in vec2 vertTexCoord0;

out vec4 fragColor;

uniform float uTime = 0.0;

void main( void )
{
    fragColor = vec4( vertTexCoord0.x, vertTexCoord0.y, 0, 1 );
}
