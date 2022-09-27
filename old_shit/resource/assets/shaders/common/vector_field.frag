#version 410

uniform sampler2D buffer;

in vec2 vertTexCoord0;

uniform vec3 iResolution;

out vec4 fragColor;

void main(void){
    vec2 C = vertTexCoord0 * iResolution.xy;
    fragColor = vec4(vec3(texelFetch(buffer,ivec2(C),0).z),1.0);
}
