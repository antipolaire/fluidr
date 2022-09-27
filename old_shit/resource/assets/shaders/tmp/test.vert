#version 410

// From our C++ code, we will tell the
// shader how many seconds have elapsed.
uniform float uTime = 0.0;

// Adjust for different aspect ratios.
uniform float uAspectRatio = 1.0;

// The GPU rasterizer has interpolated the
// texture coordinates, which now have the
// correct values for this pixel (fragment).
in vec2 vertTexCoord0;

// The output is an RGBA color, even if we
// only want grayscale.
out vec4 fragColor;

void main( void )
{
    // Start with black.
    fragColor = vec4( 1, 0, 1, 1 );
    
    // Be creative :)
    vec2 uv = vertTexCoord0 * 2.0 - 1.0;
    uv.x *= uAspectRatio;
    
    float circular = dot( uv, uv );
    float polar = atan( uv.y, uv.x );
    
    fragColor.r += 0.5 + 0.5 * cos( 4.0 * circular + 0.2 * uTime );
    fragColor.r *= 0.5 + 0.5 * cos( 5.0 * polar - 0.5 * uTime );
}
