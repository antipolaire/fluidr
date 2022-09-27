#version 410

uniform sampler2D iChannel0;

in highp vec2 vertTexCoord0;

uniform vec3 iResolution;

uniform vec2 iMouse;

out vec4 fragColor;

uniform float uTime = 0.0;

vec4 t (vec2 v, int a, int b) {return texture(iChannel0,((v+vec2(a,b))/iResolution.xy));}
vec4 t (vec2 v) {return texture(iChannel0,(v/iResolution.xy));}
float area (vec2 a, vec2 b, vec2 c) { // area formula of a triangle from edge lengths
    float A = length(b-c), B = length(c-a), C = length(a-b), s = 0.5*(A+B+C);
    return sqrt(s*(s-A)*(s-B)*(s-C));
}

void main()
{
    fragColor.w = 1.0;
    
    vec2 U = vertTexCoord0 * iResolution.xy;
    
    vec2 ur = iResolution.xy;
    if (uTime < 1||U.x < 3.||ur.x-U.x < 3.) {
        fragColor = vec4(0.1,0,0,0);
        
    } else {
        vec2 v = U,
        A = v + vec2( 1, 1),
        B = v + vec2( 1,-1),
        C = v + vec2(-1, 1),
        D = v + vec2(-1,-1);
        for (int i = 0; i < 2; i++) {
            vec2 tmp = t(v).xy;
            v -= tmp;
        }
        vec4 me = t(v);
        for (int i = 0; i < 3; i++) {
            vec2 tmp = t(v).xy;
            v -= tmp;
        }
        me.zw = t(v).zw;
        for (int i = 0; i < 9; i++) {
            A += t(A).xy;
            B += t(B).xy;
            C += t(C).xy;
            D += t(D).xy;
        }
        vec4 n = t(v,0,1),
        e = t(v,1,0),
        s = t(v,0,-1),
        w = t(v,-1,0);
        vec4 ne = .25*(n+e+s+w);
        me = mix(me,ne,vec4(0.06,0.06,1,0.0));
        me.z  = me.z + (area(A,B,C)+area(B,C,D)-4.);
        vec4 pr = vec4(e.z,w.z,n.z,s.z);
        me.xy = me.xy + vec2(pr.x-pr.y, pr.z-pr.w)/ur;
        
        U = U-vec2(0.4,0.5)*ur;
        float an = -iMouse.x/ur.x,
        co = cos(an), si = sin(an);
        U.xy = mat2(co,-si,si,co)*U.xy;
        U.x*=0.125;
        U.y += (iMouse.y/ur.y)*U.x*U.x;
        
        me.xyw *= step(6.,length(U));
        fragColor = me;
        fragColor.xyz = clamp(fragColor.xyz, -vec3(.5,.5,40.), vec3(.5,.5,40.));
    }
}
