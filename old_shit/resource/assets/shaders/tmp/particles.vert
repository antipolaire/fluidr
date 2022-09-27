#version 410

// Automatically provided by Cinder.
uniform mat4 ciModelViewProjection;

layout(location = 0) in vec4 iPositionVelocity;  // xy = position, zw = velocity

void main(void)
{
	// Just output the position and a size.
	gl_PointSize = 8.0;
	gl_Position = ciModelViewProjection * vec4( iPositionVelocity.xy, 0.0, 1.0 );
}
