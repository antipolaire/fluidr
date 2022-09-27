#version 410

precision highp float;

#include "texture_utils.glsl"

// Read: http://developer.download.nvidia.com/books/HTML/gpugems/gpugems_ch38.html
// Code: https://http.download.nvidia.com/developer/GPU_Gems/CD_Image/Index.html

uniform sampler2D iChannel0;

uniform sampler2D iChannel3;

in vec2 vertTexCoord0;

uniform vec3 iResolution;

out vec4 fragColor;

uniform float uTime = 0.0;

uniform vec2 iMouse;

//void advect(float2 coords   : WPOS,   // grid coordinates
//
//            out float4 xNew : COLOR,  // advected qty
//
//            uniform float timestep,
//            uniform
//            float rdx,        // 1 / grid scale
//
//            uniform
//            samplerRECT u,    // input velocity
//
//            uniform
//            samplerRECT x)    // qty to advect
//{
//    // follow the velocity field "back in time"
//
//    float2 pos = coords - timestep * rdx * f2texRECT(u, coords);
//
//    // interpolate and write to the output fragment
//    xNew = f4texRECTbilerp(x, pos);
//}

// Runge-Kutta 4 backward advection

#define h 2.

// #define has_obstacle true

// Nice documentation: https://www.haroldserrano.com/blog/visualizing-the-runge-kutta-method

vec2 RK4(vec2 p){
    vec2 k1 = texelFetch(iChannel0,ivec2(clamp(p, vec2(0), iResolution.xy-vec2(1))),0).xy;
    vec2 k2 = texelFetch(iChannel0,ivec2(clamp(p-0.5*h*k1, vec2(0), iResolution.xy-vec2(1))),0).xy;
    vec2 k3 = texelFetch(iChannel0,ivec2(clamp(p-0.5*h*k2, vec2(0), iResolution.xy-vec2(1))),0).xy;
    vec2 k4 = texelFetch(iChannel0,ivec2(clamp(p-h*k3, vec2(0), iResolution.xy-vec2(1))),0).xy;
    return h/3.*(0.5*k1+k2+k3+0.5*k4);
}

vec2 RK42(vec2 p){
    vec2 k1 = texelFetch(iChannel0,ivec2(p),0).xy;
    vec2 k2 = texelFetch(iChannel0,ivec2(p.x+h*.5,p.y+0.5*h*k1),0).xy;
    vec2 k3 = texelFetch(iChannel0,ivec2(p.x+h*.5,p.y+0.5*h*k2),0).xy;
    vec2 k4 = texelFetch(iChannel0,ivec2(p.x,p.y+h*k3),0).xy;
    return h/6.*(k1+2.*k2+2.*k3+k4);
}

void render(in vec2 C )
{
    vec2 resolution_xy = iResolution.xy;

    vec4 buf = texelFetch(iChannel0,ivec2(clamp(C-RK42(C),vec2(0),iResolution.xy-vec2(1))),0); // advect

    // vec4 buf = texelFetch(iChannel0,ivec2(C)-ivec2(RK4(C)),0); // advect

    vec2 velocity_xy = buf.xy;
    float density_xy = buf.z;
    
    // set boundary velocity
    if( abs(resolution_xy.x-C.x)<=2.0||
        abs(C.x)<=2.0||
        abs(resolution_xy.y-C.y)<=2.0||
        abs(C.y)<=2.0){
        velocity_xy = vec2(0,0);
    }
    
    // set density for floating "bars" from the left
    if(0.<C.x&&C.x<2.){
        if(cos(C.y*0.4)>0.6){
            density_xy = 1.0;
        }else{
            density_xy = 0.;
        }
    }

    // Read: iMouse stores drag and drop coordinates which make it
    // possible to calculate the length: https://shadertoyunofficial.wordpress.com/2016/07/20/special-shadertoy-features/

//    vec2 obstacle_coordinates_xy = vec2(100*cos(uTime), 100*sin(uTime))+iResolution.xy*0.5;
//
//    vec2 m = (obstacle_coordinates_xy-0.5*resolution_xy)/resolution_xy.y;
//    if(length(obstacle_coordinates_xy)<0.01){
//        m = vec2(-0.5,0.);
//    }
//    vec2 velocity_mouse = m-vec2(texelFetch(iChannel0,ivec2(C)+ivec2(1,0),0).a,
//                                 texelFetch(iChannel0,ivec2(C)+ivec2(0,1),0).a); // mouse velocity

    if(uTime<5){
        velocity_xy = vec2(1.0,0.0);
        //velocity_mouse = vec2(0);
    }

    /* Velocity of our dot!*/
//    if(has_obstacle && distance(obstacle_coordinates_xy.xy,C)<10){
//        velocity_xy = 100.*velocity_mouse;
//        density_xy = min(1.,length(velocity_mouse)*100.);
//    }

    //float mxy = m.x;
    fragColor = vec4(velocity_xy,density_xy, 1.0);
}

void main(void){
    render(vertTexCoord0*iResolution.xy);
}
