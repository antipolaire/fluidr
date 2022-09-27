#version 410

// We'd like to know the size of the window.
uniform vec2 uWindowSize;

// We'd like to sample the dynamic background.
uniform sampler2D uTexBackground;

// Automatically provided by Cinder.
uniform mat4 ciModelViewProjection;

// Input attributes: our position and velocity.
layout(location = 0) in vec4 iPositionVelocity;

// Output attributes: should be the same.
layout(location = 0) out vec4 oPositionVelocity;

void main( void )
{
	// Start by copying the current data.
	oPositionVelocity = iPositionVelocity;

	// Convert particle position to a normalized texture coordinate,
	// so that we can sample the background texture. This is simple 
	// in our case, because the texture has the same size as the window.
	vec2 coord = oPositionVelocity.xy / uWindowSize;

	// Sample the background image four times and
	// determine the slope. 
	float top = textureOffset( uTexBackground, coord, ivec2(0, -1) ).r;
	float left = textureOffset( uTexBackground, coord, ivec2(-1, 0) ).r;
	float bottom = textureOffset( uTexBackground, coord, ivec2(0, 1) ).r;
	float right = textureOffset( uTexBackground, coord, ivec2(1, 0) ).r;
	vec2 slope = vec2( right - left, bottom - top );

	// Update velocity. Particle will slow down when going uphill.
	oPositionVelocity.zw = 0.998 * oPositionVelocity.zw - 1.0 * slope;

	// Update position.
	oPositionVelocity.xy += oPositionVelocity.zw;

	// Make sure the particle does not leave the window.
	if( oPositionVelocity.x < 0 )
		oPositionVelocity.x += uWindowSize.x;
	else if( oPositionVelocity.x >= uWindowSize.x )
		oPositionVelocity.x -= uWindowSize.x;

	if( oPositionVelocity.y < 0 )
		oPositionVelocity.y += uWindowSize.y;
	else if( oPositionVelocity.y >= uWindowSize.y )
		oPositionVelocity.y -= uWindowSize.y;
}