#version 410

precision highp float;

#include "texture_utils.glsl"

uniform sampler2D iChannel0;

in vec2 vertTexCoord0;

uniform vec3 iResolution;

out vec4 fragColor;

uniform float uTime = 0.0;

// Compute divergence

void render( in vec2 C ){
    float wLx = texelFetch(iChannel0,ivec2(C)+ivec2(-1,0),0).x; // texelFetch(iChannel0,clamp(ivec2(C-vec2( 1,0)), ivec2(0), ivec2(iResolution.xy)-1),0).x;
    float wRx = texelFetch(iChannel0,ivec2(C)+ivec2(1,0),0).x; // texelFetch(iChannel0,clamp(ivec2(C+vec2( 1,0)), ivec2(0), ivec2(iResolution.xy)-1),0).x;
    float wBy = texelFetch(iChannel0,ivec2(C)+ivec2(0,-1),0).y; // texelFetch(iChannel0,clamp(ivec2(C-vec2( 0,1)), ivec2(0), ivec2(iResolution.xy)-1),0).y;
    float wTy = texelFetch(iChannel0,ivec2(C)+ivec2(0,1),0).y; // texelFetch(iChannel0,clamp(ivec2(C+vec2( 0,1)), ivec2(0), ivec2(iResolution.xy)-1),0).y;

    float div = 0.1 * ((wRx-wLx)+(wTy-wBy));
    fragColor = vec4(div,0,0,1);
}

//void divergence(half2 coords  : WPOS,   // grid coordinates
//
//                out
//                half4 div : COLOR,  // divergence
//
//                uniform half halfrdx,   // 0.5 / gridscale
//                
//                uniform samplerRECT w)  // vector field
//{
//    half4 wL = h4texRECT(w, coords - half2(1, 0));
//    half4 wR = h4texRECT(w, coords + half2(1, 0));
//    half4 wB = h4texRECT(w, coords - half2(0, 1));
//    half4 wT = h4texRECT(w, coords + half2(0, 1));
//
//    div = halfrdx * ((wR.x - wL.x) + (wT.y - wB.y));
//}

void main(void){
    render(vertTexCoord0*iResolution.xy);
}
