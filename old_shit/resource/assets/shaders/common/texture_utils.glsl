vec4 textureMS(in sampler2D tex, in vec2 texCoord, in vec2 pixelOffset)
{
  ivec2 textureSize = textureSize(tex,0);
  ivec2 texelCoords = ivec2(textureSize * clamp(texCoord, 0.0, 1.0) + pixelOffset);
  return texelFetch(tex, texelCoords, 0);
}