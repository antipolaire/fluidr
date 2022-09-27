#version 410


uniform sampler2D iChannel0;

uniform sampler2D iChannel1;

uniform sampler2D iChannel2;

uniform sampler2D iChannel3;

in vec2 vertTexCoord0;

uniform vec3 iResolution;

uniform vec2 iMouse;

out vec4 fragColor;

uniform float uTime = 0.0;

uniform int renderChannel = 0;

void main(void){
    vec2 C = vertTexCoord0*iResolution.xy;
    if(distance(iMouse.xy,C)<10){
        fragColor =  vec4(0.3,0.5,0.7,1);
    }else{
        switch(renderChannel){
            case 0:
                fragColor = vec4(texture(iChannel0,vertTexCoord0).z);
                break;
            case 1:
                fragColor = vec4(texture(iChannel1,vertTexCoord0).z);
                break;
            case 2:
                fragColor = vec4(texture(iChannel2,vertTexCoord0).z);
                break;
            case 3:
                fragColor = vec4(texture(iChannel3,vertTexCoord0).z);
                break;
            default:
                fragColor = vec4(0.0,0.0,0.0,1.0);
        }
    }
}
