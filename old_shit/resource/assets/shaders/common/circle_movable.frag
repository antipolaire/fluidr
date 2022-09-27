#version 410

uniform sampler2D iChannel3;

in vec2 vertTexCoord0;

uniform vec3 iResolution;

uniform vec2 _iMouse;

out vec4 fragColor;

uniform float uTime = 0.0;

void main(void){
    vec2 C = vertTexCoord0*iResolution.xy;

    vec2 iMouse = vec2(100*cos(uTime), 100*sin(uTime))+iResolution.xy*0.5;

    if(distance(iMouse.xy,C)<10){
        fragColor =  vec4(0.3,0.5,0.7,1);
    }else{
        fragColor = vec4(texelFetch(iChannel3,clamp(ivec2(C), ivec2(0), ivec2(iResolution.xy)-1),0).z);
    }
}
