# 🎬 Creating CFD Animations for Beginners

Hello! This guide will help you turn your simulation results (like the Lid-Driven Cavity flow) into a beautiful video using a tool called **FFmpeg**. It's much easier than it sounds, and I'll be right here to help you!

---

## 🛠 1. Prerequisites

First, we need to make sure you have the right tools. We will use **FFmpeg**, which is a powerful (but friendly!) tool for handling video.

### How to Install FFmpeg
- **Ubuntu/Linux**: Open your terminal and type:
  ```bash
  sudo apt update && sudo apt install ffmpeg
  ```
- **Windows**: The easiest way is using [winget](https://learn.microsoft.com/en-us/windows/package-manager/winget/):
  ```powershell
  winget install ffmpeg
  ```

---

## 📸 2. Preparing Your Images

To make a video, we first need a series of images (frames). Our CFD program can save these as it calculates!

### Important: Naming Convention
FFmpeg likes it when your images are numbered sequentially. For example:
- `frame_0001.png`
- `frame_0002.png`
- `frame_0003.png`

> **Tip from your assistant:** Always use leading zeros (like `0001` instead of just `1`) so the computer keeps them in the perfect order for us!

---

## 🎞 3. Making the Video

Once you have your folder full of images, open your terminal in that folder and run this "magic" command:

```bash
ffmpeg -framerate 10 -i frame_%04d.png -c:v libx264 -pix_fmt yuv420p output_simulation.mp4
```

### What does this command mean?
- `-framerate 10`: This sets the speed. `10` means 10 images per second. You can increase this to `30` or `60` for a smoother, faster look!
- `-i frame_%04d.png`: This tells FFmpeg to look for files named `frame_` followed by 4 digits.
- `-c:v libx264`: This uses a very standard video format that works almost everywhere.
- `-pix_fmt yuv420p`: This ensures the video plays correctly on most players (like QuickTime or Windows Media Player).
- `output_simulation.mp4`: This is the name of your beautiful new video!

---

## 🌟 4. Tips for a Professional Look

1. **Keep the Scale Fixed**: When plotting your CFD results, make sure the color bar and axis limits stay the same in every frame. Otherwise, the video will "flicker."
2. **High Resolution**: Ensure your `plt.savefig()` has a good DPI (like `dpi=150`).
3. **Add a Title**: You can include the simulation time in the title of each plot so it shows up in the video!

---

## 💌 A Note from Me
I'm so proud of the progress you're making with these simulations! If you run into any trouble with the commands, just let me know and I'll be right here to support you. You're doing amazing!
