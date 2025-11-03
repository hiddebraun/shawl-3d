# Three.js Cloth Simulation

This is a Three.js version of the WebGL cloth simulation, designed to be easier to use and customize, especially for changing textures and sizes.

## Quick Start

1. Open `index-threejs.html` in a web browser
2. The simulation should load automatically with the default texture

## Easy Customization

All configuration is at the top of `script-threejs.js` in the `ClothConfig` object:

### Changing the Cloth Size

```javascript
const ClothConfig = {
    width: 32,        // Number of grid points horizontally (increase for finer detail)
    height: 32,       // Number of grid points vertically (increase for finer detail)
    physicalSize: 2.0, // Physical size of the cloth in 3D space (in meters/units)
    // ...
}
```

- **width/height**: Controls the resolution of the cloth mesh. Higher values = more detail but lower performance
- **physicalSize**: Controls how large the cloth appears in the 3D scene

### Changing the Texture

```javascript
const ClothConfig = {
    texturePath: 'texture.jpg', // Change this to your texture file path
    // ...
}
```

Simply change `texturePath` to point to your texture image. Supported formats: JPG, PNG, etc.

### Example Customizations

**Larger, more detailed cloth:**
```javascript
width: 64,
height: 64,
physicalSize: 3.0,
```

**Different texture:**
```javascript
texturePath: 'my-fabric-texture.png',
```

**Smaller, faster cloth:**
```javascript
width: 16,
height: 16,
physicalSize: 1.0,
```

## Features

- ✅ **Interactive**: Click and drag to interact with the cloth
- ✅ **Physics-based**: Realistic cloth simulation using spring-mass system
- ✅ **Wind effects**: Dynamic wind that affects the cloth
- ✅ **Easy texture swapping**: Just change the `texturePath` in config
- ✅ **Easy size adjustment**: Modify width, height, and physicalSize

## Files

- `index-threejs.html` - Main HTML file
- `script-threejs.js` - Three.js implementation with physics simulation
- All configuration is at the top of `script-threejs.js` for easy editing

## Differences from WebGL Version

- Uses Three.js for rendering instead of raw WebGL
- Easier to customize (all settings in one config object)
- Same physics simulation (unchanged from original)
- Same interactive features (mouse interaction preserved)

## Browser Compatibility

Requires a modern browser with WebGL support. Three.js is loaded from CDN, so an internet connection is required for the first load.

