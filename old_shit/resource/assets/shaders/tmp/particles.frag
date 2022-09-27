#version 410

// The output is an RGBA color.
out vec4 fragColor;

void main(void)
{
	// Make the particles green.
	fragColor.rgb = vec3( 0, 0.25, 0.1 );

	// With a little trick, we can make them round.
	vec2 uv = gl_PointCoord.xy * 2.0 - 1.0;
	float r = dot( uv, uv );
	fragColor.a = 1.0 - smoothstep( 0.45, 0.55, r );
}