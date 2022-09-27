#version 410

uniform sampler2D iChannel0;

in vec2 vertTexCoord0;

out vec4 fragColor;

uniform float uTime = 0.0;

void main( void )
{
    vec3 texColor = texture(iChannel0,vertTexCoord0).rgb;
    fragColor = vec4(texColor.x,texColor.y,texColor.z,1.0);
    
    vec2 uv = vertTexCoord0 * 2.0 - 1.0;
    float circular = dot( uv, uv );
    fragColor.b += 0.5 + 0.5 * cos( 40.0 * circular + 2.0 * uTime );
}
