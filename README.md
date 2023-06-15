# rela-tracer: A Special Relativistic Raytracer

rela-tracer is a tool for rendering objects that move at speeds close to the speed of light. It uses stochastic, physics-based lighting to produce realistic images.

## Table of Contents <!-- omit in toc -->

- [How it Works](#how-it-works)
- [Example scenes](#example-scenes)
  - [Glowing Spheres](#glowing-spheres)
  - [Headlight Absorption](#headlight-absorption)
  - [Dielectric Delays](#dielectric-delays)
  - [Box Array](#box-array)
  - [Boxes in Mirror](#boxes-in-mirror)
  - [Moving Box](#moving-box)
  - [Metal Balls](#metal-balls)
  - [Mirror Box](#mirror-box)
- [Installation](#installation)
- [Usage](#usage)
- [Limitations](#limitations)
- [License](#license)

## How it Works

Each scene in rela-tracer consists of a collection of objects and camera data. Every object has its own inertial reference frame, in which it is completely stationary.

To render an image, rela-tracer follows a similar process to other raytracing programs. It calculates the paths of light rays entering the camera for each pixel in the image. Then, it follows these rays backward in time to determine the closest intersection with an object in the scene.
Once an object is hit rela-tracer adds any light emitted by the object to the pixel color (taking into account special relativistic effects), and performs a _stochastic inverse scatter_: It randomly picks a ray which could have scattered off in the direction of the ray in question, then traces said new ray. This process is terminated after a specified number of bounces ("depth").

Rela-tracer implements the physics of special relativity by breaking up the computation of intersections between a light ray and an object into three steps. First, it transforms the light ray into the restframe of the object. Then, it computes intersection data for the light ray and object in this frame. Finally, it transforms the intersection point and time back to the standard frame.

Because of the relativistic Doppler effect, light rays not only change direction but also change color and brightness when transformed from one frame to another. Rela-tracer takes this into account by keeping track of the cumulative relativistic transform to be applied to light coming from hit objects. This also includes absorption effects.

## Example scenes

The following example scenes are available in this project:

### Glowing Spheres

Three identical spheres emitting green light are placed in front of a mirror. The top sphere moves towards the camera, the bottom sphere away from it. As a result, the top sphere appears blue-shifted and the bottom sphere red-shifted. In the mirror, the effect is reversed. The mirror is tinted in order for it to be visible against the background.

https://github.com/c-weis/rela-tracer/assets/34036773/f18a7491-4ff9-4a58-84e6-5842361295d2

### Headlight Absorption
Three diffusely reflecting spheres are placed inside a box. The right/left panels emit red/green light, respectively. The top sphere moves towards the right, while the bottom sphere moves towards the left. This scene showcases the headlight effect for absorption: just as a fast moving light emitter emits most of its light along its direction of movement, so do fast moving objects receive most of their light along that direction.

https://github.com/c-weis/rela-tracer/assets/34036773/53e3649b-7d6f-4b90-9786-c19e0d99f0f6

### Dielectric Delays

Two light-emitting balls travel downwards. The view towards the left ball is obstructed by two boxes made from dielectric materials: the top one is made from glass, while the bottom one is made from diamond. Light travels slower in dielectric media, and this scene showcases the resulting delay. Further, the Doppler effect is visible again both across the range of movement of the balls and in the internal reflections of the balls in sides of the boxes.

https://github.com/c-weis/rela-tracer/assets/34036773/eb9e22d8-e2ad-4d22-98da-7ff661f560e9

### Box Array

A vertical array of boxes is moving horizontally at varying speeds. They are illuminated by ambient white light. Their material is set to absorb short wavelengths, making the resting box appear yellow. The scene showcases the Doppler effect on absorption properties.

https://github.com/c-weis/rela-tracer/assets/34036773/5427b28f-e16b-4393-98c6-1949f2c66eec

### Boxes in Mirror

A row of boxes move vertically at varying speeds above of a mirror. The fast-moving boxes appear distorted, and the mirror reflects their images with an apparent time delay due to longer light travel time.

https://github.com/c-weis/rela-tracer/assets/34036773/b5e5f4b1-1ec2-4f07-9f9f-8b0e4d012dd5

### Metal Balls

A static scene featuring three reflective balls, two of which are tinted. The balls are placed in an environment with a ground plane and a skylight. The image demonstrates reflections and the interaction of the two different absorption spectra.

![metal_balls](https://github.com/c-weis/rela-tracer/assets/34036773/f6b5b263-caed-430d-adab-6c005c9feb81)

### Mirror Box

A metallic sphere and a small spherical light source are placed inside a box with reflective sides. The right panel of the box is tinted red, while the top emits light.

![mirror_box](https://github.com/c-weis/rela-tracer/assets/34036773/e6277ada-2921-458c-b656-9885ce864319)

## Installation

Clone the repository and run make. This creates the executable _rela-tracer_.
E.g.

    ```bash
    git clone https://github.com/c-weis/rela-tracer.git
    cd rela-tracer
    make
    ```

This should create an executable _rela-tracer_.

## Usage

When called without command line arguments,
_rela-tracer_ renders the scene [Headlight Absorption](#headlight-absorption) and outputs to the file `output.bmp`.
Single image renders are output to the `images` subfolder of the runtime directory, while frames of film renders are output to the `vidframes` subfolder.
The example scenes may be rendered by calling:

```bash
./rela-tracer --scene_name <string>
```

where _\<string\>_ is one of

-   "headlight_absorption": [Headlight Absorption](#headlight-absorption)
-   "glowing_spheres": [Glowing Spheres](#glowing-spheres)
-   "dielectric_delays": [Dielectric Delays](#dielectric-delays)
-   "box_array": [Box Array](#box-array)
-   "boxes_in_mirror": [Boxes in Mirror](#boxes-in-mirror)
-   "moving_box": [Moving Box](#moving-box)
-   "metal_balls": [Metal Balls](#metal-balls)
-   "mirror_box": [Mirror Box](#mirror-box)
-   "custom_scene": empty scene, to be customised in code

You can modify the scenes or create new ones by editing the `main.cc` file. The `custom_scene` function is an empty template to start from.

In addition, _rela-tracer_ accepts command line arguments that allow rendering scenes with different settings. Here is the complete list of available arguments:

```
-sn, --scene_name <string>
    Name of the example scene (default: "headlight_absorption").

-fn, --filename <string>
    Name of the output filename (default: "output").

-t, --test
    Run tests instead of rendering (default: false).

-w, --width <int>
    Image width in pixels (default: 500).

-h, --height <int>
    Image height in pixels (default: 500).

-rpp, --rays_per_pixel <int>
    Number of rays to average over per pixel (default: 100).

-d, --depth <int>
    Maximum number of bounces to trace rays (default: 8).

-ipu, --iterations_per_update <int>
    Number of passes between image updates (default: 1).

-f, --film
    Render frames of a film (default: false).

-st, --start_time <float>
    Time of the picture in standard frames (default: -1.0).

-et, --end_time <float>
    Time of the last frame for a film in standard frames (default: 3.0).

-dt, --d_time <float>
    Time between frames for a film in standard frames (default: 0.1).

-po, --preview_only
    Generate only preview picture(s) (default: false).

-sf, --start_frame <int>
    First frame to render (default: 0).

-ef, --end_frame <int>
    Last frame to render, -1 = all remaining (default: -1).

-ci, --camera_index <int>
    Index of the camera, a scene may contain several (default: 0).

-rf, --rescale_factor <float>
    Multiplier for RGB values, -1 = estimate it (default: -1).

```

Here are some example use cases:

1. Render [Headlight absorption](#headlight-absorption) with a custom output filename:

    ```bash
    rela-tracer --scene_name my_scene --filename custom_output
    ```

2. Render a preview picture of a scene named "example_scene" with a width of 800 pixels and 50 rays per pixel:

    ```bash
    rela-tracer -sn example_scene -w 800 -rpp 50 -po
    ```

3. Run tests instead of rendering a scene named "test_scene":

    ```bash
    rela-tracer -t --scene_name test_scene
    ```

4. Render frames of a film named "my_film" with a start time of 0.5, an end time of 5.0, and a time increment of 0.2 between frames:

    ```bash
    rela-tracer -f -sn my_film -st 0.5 -et 5.0 -dt 0.2
    ```

5. Render a scene named "advanced_scene" with a custom rescale factor of 2.5:
    ```bash
    rela-tracer --scene_name advanced_scene --rescale_factor 2.5
    ```

These examples demonstrate the usage of different command line arguments to customize the rendering process according to specific requirements. Remember to use either the short form (single dash) or the long form (double dash) consistently for each argument.

```bash
rela-tracer --width 800 --height 600 --depth 12 --filename headlights --film
```

This command will render a film with a width of 800 pixels, a height of 600 pixels, a maximum ray depth of 12, and save the output frames with the prefix "scene_output".

## Limitations

By its nature, rela-tracer cannot render rotating objects because a full treatment of rotating objects requires the framework of general relativity, which involves non-Euclidean geometry.

The stochastic character of the simulation often leads to grainy images, especially for

-   low number of rays per pixel
-   diffuse or fuzzy surfaces
-   sparsely lit scenes.

Further, we make a number of approximations which may impact the resulting images: the conversion from color spectra in terms of wavelengths into RGB-values proceeds via multimodal Gaussian fits to the CIE 1931 XYZ color functions. The error introduced is on the order of 1%.

## License

The ray tracing project is licensed under the MIT License. You can find the full license text in the `LICENSE` file.
