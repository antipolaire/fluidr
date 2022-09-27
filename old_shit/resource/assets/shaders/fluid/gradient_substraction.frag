#version 410

precision highp float;

#include "texture_utils.glsl"

uniform sampler2D iChannel2;

uniform sampler2D iChannel0;

in vec2 vertTexCoord0;

uniform vec3 iResolution;

out vec4 fragColor;

uniform float uTime = 0.0;

// Gradient subtraction

void render( in vec2 C ){

    // iChannel2 = pressure
    float pL = texelFetch(iChannel2,ivec2(C)+ivec2(-1,0),0).x;
    float pR = texelFetch(iChannel2,ivec2(C)+ivec2(1,0),0).x;
    float pB = texelFetch(iChannel2,ivec2(C)+ivec2(0,-1),0).x;
    float pT = texelFetch(iChannel2,ivec2(C)+ivec2(0,1),0).x;

    // iChannel0 = velocity
    vec4 uNew = texelFetch(iChannel0,ivec2(C),0);

    vec2 grad = 0.1 * vec2(pR-pL,pT-pB);
    fragColor = vec4(uNew.xy-grad,uNew.z,1);
}

//void gradient(half2 coords   : WPOS,   // grid coordinates
//
//              out half4 uNew : COLOR,  // new velocity
//
//              uniform half halfrdx,    // 0.5 / gridscale
//
//              uniform samplerRECT p,   // pressure
//
//              uniform samplerRECT w)   // velocity
//{
//    half pL = h1texRECT(p, coords - half2(1, 0));
//    half pR = h1texRECT(p, coords + half2(1, 0));
//    half pB = h1texRECT(p, coords - half2(0, 1));
//    half pT = h1texRECT(p, coords + half2(0, 1));
//
//    uNew = h4texRECT(w, coords);
//    uNew.xy -= halfrdx * half2(pR - pL, pT - pB);
//}

void main(void){
    render(vertTexCoord0*iResolution.xy);
}
