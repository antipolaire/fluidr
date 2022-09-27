#version 410

// Automatically provided by Cinder.
uniform mat4 ciModelViewProjection;

// Input attributes automatically provided by Cinder.
in vec4 ciPosition;
in vec2 ciTexCoord0;

// Output attributes which we will have to provide ourselves.
out vec2 vertTexCoord0;

void main( void )
{
    // Output the vertex coordinate. The GPU will
    // interpolate it for us, so that it has the
    // correct value for every pixel.
    vertTexCoord0 = ciTexCoord0;
    
    // OpenGL requires us to transform every vertex
    // to so-called normalized device coordinates.
    // This sounds harder than it is:
    gl_Position = ciModelViewProjection * ciPosition;
}
