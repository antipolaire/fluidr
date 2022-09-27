#version 410

uniform sampler2D iChannel0;

in highp vec2 vertTexCoord0;

uniform vec3 iResolution;

uniform vec2 iMouse;

out vec4 fragColor;

uniform float uTime = 0.0;

void main()
{
    
    fragColor.a = 1.0;
    
    vec2 U = vertTexCoord0 * iResolution.xy;
    
    vec4 g = texture(iChannel0,U/iResolution.xy);
 
    fragColor.rgb = vec3(g.a);
    
    U = U-vec2(0.4,0.5)*iResolution.xy;
    float an = -iMouse.x/iResolution.x,
    co = cos(an), si = sin(an);
    U.xy = mat2(co,-si,si,co)*U.xy;
    U.x*=0.125;
    U.y += (iMouse.y/iResolution.y)*U.x*U.x;
    if (length(U)<6.) {
        fragColor.xyz = vec3(sin(uTime)*0.5+0.5,0.5,1.);
    }
}
