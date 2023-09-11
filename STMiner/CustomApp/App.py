import tkinter as tk
from tkinter import filedialog, messagebox
from tkinter import simpledialog

from PIL import Image, ImageTk, ImageDraw
from PIL import ImageFile

# Allow loading large images
Image.MAX_IMAGE_PIXELS = None
ImageFile.LOAD_TRUNCATED_IMAGES = True


class App:
    def __init__(self):
        self.img = None
        self.img_on_canvas = None

        # Root windows
        self.root = tk.Tk()
        self.root.title("Region Annotation")

        # Create Canvas
        self.canvas = tk.Canvas(self.root, width=800, height=800)
        self.canvas.pack()
        self.canvas.bind("<Button-1>", self._annotate)
        self.canvas.bind("<B1-Motion>", self._continue_drawing)
        self.canvas.bind("<ButtonRelease-1>", self._stop_drawing)
        self.canvas.bind("<Button-3>", self._annotate)

        # Create menu
        self.menu = tk.Menu(self.root)
        self.root.config(menu=self.menu)
        self.file_menu = tk.Menu(self.menu)
        self.menu.add_cascade(label="File", menu=self.file_menu)
        self.file_menu.add_command(label="Load image", command=self.open_image)
        self.file_menu.add_command(label="Save image", command=self.save_image)

        self.op_menu = tk.Menu(self.menu)
        self.menu.add_cascade(label="Edit", menu=self.op_menu)
        self.op_menu.add_command(label="Brush +", command=self._upper)
        self.op_menu.add_command(label="Brush -", command=self._lower)
        self.op_menu.add_command(label="Cut image", command=self._cut)
        self.op_menu.add_command(label="Clean up annotations", command=self.clear_annotations)
        self.op_menu.add_command(label="Reset", command=self.reset)

        self.help_menu = tk.Menu(self.menu)
        self.menu.add_cascade(label="Help", menu=self.help_menu)
        self.help_menu.add_command(label="Help", command=self.show_help)
        self.help_menu.add_command(label="About", command=self.show_about)

        self.drawing = False
        self.radius = 5
        self.left = None
        self.right = None
        self.top = None
        self.bottom = None
        self.origin_img = None

    def _cut(self):
        user_input = simpledialog.askstring("Input edge",
                                            "Input 'left,right,top,bottom' (Comma separated):")
        if user_input:
            params = user_input.split(',')
            if len(params) == 4:
                self.left, self.right, self.top, self.bottom = params
                self.root.title("Cutting image... Please wait.")
                self.img = self.origin_img.crop((int(self.left),
                                                 int(self.top),
                                                 int(self.right),
                                                 int(self.bottom)))
                self.img = self.img.resize((800, 800), Image.ANTIALIAS)
                self.img_on_canvas = ImageTk.PhotoImage(self.img)
                # self.canvas.delete("all")
                self.canvas.create_image(0, 0,
                                         anchor=tk.NW,
                                         image=self.img_on_canvas)
                self.root.title("Region Annotation")
            else:
                messagebox.showinfo("Error", 'Format Error!')

    def _upper(self):
        self.radius += 1

    def _lower(self):
        if self.radius > 1:
            self.radius -= 1

    def open_image(self):
        file_path = filedialog.askopenfilename()
        if file_path:
            self.root.title("Loading image... Please wait.")
            self.origin_img = Image.open(file_path)
            self.img = self.origin_img.resize((800, 800), Image.ANTIALIAS)
            self.img_on_canvas = ImageTk.PhotoImage(self.img)
            self.canvas.create_image(0, 0, anchor=tk.NW, image=self.img_on_canvas)
            self.root.title("Region Annotation")

    def _draw(self, event):
        x, y = event.x, event.y
        draw = ImageDraw.Draw(self.img)
        draw.ellipse((x - self.radius, y - self.radius, x + self.radius, y + self.radius),
                     fill="red",
                     outline="red")
        self.canvas.create_oval(x - self.radius, y - self.radius, x + self.radius, y + self.radius,
                                fill="red",
                                outline="red")

    def _annotate(self, event):
        self.drawing = True
        self._draw(event)

    def _continue_drawing(self, event):
        if self.drawing:
            self._draw(event)

    def _stop_drawing(self, event):
        self.drawing = False

    def clear_annotations(self):
        self.canvas.delete("all")
        self.canvas.create_image(0, 0,
                                 anchor=tk.NW,
                                 image=self.img_on_canvas)

    def reset(self):
        self.canvas.delete("all")

    def save_image(self):
        file_path = filedialog.asksaveasfilename(defaultextension=".png",
                                                 filetypes=[('*', '.png')])
        if file_path:
            self.img.save(file_path)

    def show_help(self):
        help_text = "Usage: \n\n" \
                    "1. Use “File” -> “Load image” to load image.\n" \
                    "2. Click mouse over image to _annotate image.\n" \
                    "3. Use “File” -> “Save image” to save the annotated image.\n\n" \
                    "Check the documentation for more details."
        messagebox.showinfo("Help", help_text)

    def show_about(self):
        about_text = "Specifies the region on tissue.\n" \
                     "Author: Peisen Sun\n" \
                     "E-mail: sunpeisen@stu.xjtu.edu.cn\n" \
                     "Check the documentation for more details."
        messagebox.showinfo("About", about_text)

    def run(self):
        self.root.mainloop()
