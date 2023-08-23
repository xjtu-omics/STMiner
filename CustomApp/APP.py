import tkinter as tk
from tkinter import filedialog, messagebox
from PIL import Image, ImageTk, ImageDraw


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
        self.canvas.bind("<Button-1>", self.annotate)

        # Create menu
        self.menu = tk.Menu(self.root)
        self.root.config(menu=self.menu)
        self.file_menu = tk.Menu(self.menu)
        self.menu.add_cascade(label="File", menu=self.file_menu)
        self.file_menu.add_command(label="Load image", command=self.open_image)
        self.file_menu.add_command(label="Save image", command=self.save_image)

        self.op_menu = tk.Menu(self.menu)
        self.menu.add_cascade(label="Operation", menu=self.op_menu)
        self.op_menu.add_command(label="Clean up annotations", command=self.clear_annotations)
        self.op_menu.add_command(label="Reset", command=self.reset)

        help_menu = tk.Menu(self.menu)
        self.menu.add_cascade(label="Help", menu=help_menu)
        help_menu.add_command(label="Help", command=self.show_help)
        help_menu.add_command(label="About", command=self.show_help)

        self.help_menu = tk.Menu(self.menu)

    def open_image(self):
        file_path = filedialog.askopenfilename()
        if file_path:
            self.root.title("Loading image... please wait.")
            self.img = Image.open(file_path)
            self.img = self.img.resize((800, 800), Image.ANTIALIAS)  # adjust img size
            self.img_on_canvas = ImageTk.PhotoImage(self.img)
            self.canvas.create_image(0, 0, anchor=tk.NW, image=self.img_on_canvas)
            self.root.title("Region Annotation")

    def annotate(self, event):
        x, y = event.x, event.y
        draw = ImageDraw.Draw(self.img)
        draw.ellipse((x - 5, y - 5, x + 5, y + 5), fill="red", outline="red")
        self.canvas.create_oval(x - 5, y - 5, x + 5, y + 5, fill="red", outline="red")

    def clear_annotations(self):
        self.canvas.delete("all")
        self.canvas.create_image(0, 0, anchor=tk.NW, image=self.img_on_canvas)

    def reset(self):
        self.canvas.delete("all")

    def save_image(self):
        file_path = filedialog.asksaveasfilename(defaultextension=".png", filetypes=[('*', '.png')])
        if file_path:
            self.img.save(file_path)

    def show_help(self):
        help_text = "Usage: \n\n" \
                    "1. Use “File” -> “Load image” to load image.\n" \
                    "2. Click mouse over image to annotate image.\n" \
                    "3. Use “File” -> “Save image” to save the annotated image.\n\n" \
                    "Check the documentation for more details."
        messagebox.showinfo("Help", help_text)

    def run(self):
        self.root.mainloop()
