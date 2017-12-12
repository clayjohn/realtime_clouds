# Realtime Clouds
Experiment with generating clouds in real time on low end computer

## About
This is an experiment in recreating procedural clouds that are similar in quality to the ones in Horizon Zero Dawn

The core algorithm is very similar to what is described in the various articles on the clouds in Horizen Zero Dawn

As of the writing of this the program runs between ~15-50ms/frame on my chromebook (2.4-GHz Intel Bay Trail-M dual-core N2830)

And on my desktop computer it runs ~0.5-2ms/frame (Nvidia GTX 1050)

This project is far from finished, it will still undergo various optimizations and tweaks as well as added functionality

The goal is to make this into something that can be used in other projects

### Screenshots!
![dec7-1](https://user-images.githubusercontent.com/16521339/33801879-362f844a-dd1e-11e7-87fb-5b13e2aefc24.png)
![dec7-2](https://user-images.githubusercontent.com/16521339/33801880-38cf507c-dd1e-11e7-9d5c-ef3ca9967a9f.png)
![dec7-3](https://user-images.githubusercontent.com/16521339/33801881-3a97b8b8-dd1e-11e7-9297-3b609dc0d42d.png)


## Resources
I used many different resources while working on this, and I will use many more as it improves

#### Horizon Zero Dawn Related
[The Real Time Volumetric Cloudscapes of Horizon Zero Dawn](https://www.guerrilla-games.com/read/the-real-time-volumetric-cloudscapes-of-horizon-zero-dawn)

[Gpu Pro 7 Article](https://www.crcpress.com/GPU-Pro-7-Advanced-Rendering-Techniques/Engel/p/book/9781498742535)

[Nubis: Authoring Real-Time Volumetric Cloudscapes with the Decima Engine](https://www.guerrilla-games.com/read/nubis-authoring-real-time-volumetric-cloudscapes-with-the-decima-engine)


#### Atmospheric Scattering
Shadertoy - [Atmosphere System Test](https://www.shadertoy.com/view/XtBXDz) - valentingalea

#### General
Production volume rendering: [Book](https://www.amazon.ca/Production-Rendering-Implementation-Magnus-Wrenninge/dp/156881724X),  [2011 Course](http://magnuswrenninge.com/productionvolumerendering),  [2017 Course](https://graphics.pixar.com/library/ProductionVolumeRendering/)

[Physically Based Sky, Atmosphere and Cloud Rendering in Frostbite](https://www.ea.com/frostbite/news/physically-based-sky-atmosphere-and-cloud-rendering)
