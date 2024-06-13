import tkinter as tk
from PIL import Image, ImageTk
import math


# Material class so user can choose
class Material:
    def __init__(self, name, thermal_conductivity, specific_heat_capacity, atomic_spacing, particle_size, density,
                 image):
        self.name = name
        self.thermal_conductivity = thermal_conductivity  # W/m.K
        self.specific_heat_capacity = specific_heat_capacity  # J/kg.K
        self.atomic_spacing = atomic_spacing  # Meters
        self.particle_size = particle_size  # Also meters
        self.density = density  # Kg/m3
        self.image = image


# Particle class
class Particle:
    def __init__(self, canvas, x, y, size):
        self.canvas = canvas
        self.size = size
        self.x = x
        self.y = y
        # Initial temperature is 0
        self.temperature = 0
        self.id = canvas.create_oval(x, y, x + size, y + size, fill="#000000")

    # Update color based on temperature
    def update_color(self):
        r = min(255, max(0, int(self.temperature)))
        color = "#{:02x}0000".format(r)
        self.canvas.itemconfig(self.id, fill=color)

    def move(self, dx, dy):
        self.canvas.move(self.id, dx, dy)
        self.x += dx
        self.y += dy


class HeatTransferSimulation:
    def __init__(self, root, canvas_size, grid_size, material):
        # Canvas size is the screen, grid size is no. of particles per row / column
        self.canvas_size = canvas_size
        self.grid_size = grid_size

        # Automatically fits in particles into screen. Note that this particle.size is in pixels
        self.particle_size = (canvas_size - 50) / ((grid_size - 1) * 3 + 1)

        # Array for particles
        self.particles = []

        # This is to see if user is dragging heat source
        self.is_dragging = False

        self.material = material
        self.atomic_spacing = self.material.atomic_spacing
        self.density = self.material.density

        self.canvas = tk.Canvas(root, width=canvas_size, height=canvas_size, bg="white")
        self.canvas.pack()

        self.load_background_image()
        self.heat_source = self.canvas.create_oval(195, 195, 205, 205, fill="red")
        self.thermal_conductivity = self.material.thermal_conductivity
        self.specific_heat_capacity = self.material.specific_heat_capacity
        # Trying with 0 initial temperature for particles and 60 degree heat source
        self.initial_temperature = 0
        self.heat_source_temp = 60

        # Calculate real distance per pixel. Should be meters, not nm
        self.nm_per_pixel = self.calculate_nm_per_pixel()

        # Create grid
        self.create_particles_and_lines()

        # Initialize particle temperatures
        self.initialize_particle_temperatures()

        # Makes heat source in front of particles
        self.canvas.tag_raise(self.heat_source)

        # Buttons
        self.canvas.bind("<ButtonRelease-1>", self.stop_dragging)
        self.canvas.bind("<Motion>", self.move_heat_source)
        self.canvas.bind("<ButtonPress-1>", self.canvas_click)

    # Divides total real length of all particles + spacings / pixels used
    def calculate_nm_per_pixel(self):
        total_size_nm = self.grid_size * self.material.particle_size + (
                self.grid_size - 1) * self.material.atomic_spacing
        return total_size_nm / (self.canvas_size - 50)

    def initialize_particle_temperatures(self):
        for row in self.particles:
            for particle in row:
                particle.temperature = self.initial_temperature
                particle.update_color()

        # Apply heat source to particles

    def apply_heat_source(self):
        # Get heat source coordinates, find center of heat source
        heat_source_coords = self.canvas.coords(self.heat_source)
        heat_source_x = (heat_source_coords[0] + heat_source_coords[2]) / 2
        heat_source_y = (heat_source_coords[1] + heat_source_coords[3]) / 2

        radius_nm = 1e-7  # Radius in nanometers for applying heat

        # Calculate power of the heat source using Stefan-Boltzmann Law
        heat_source_temp_kelvin = self.heat_source_temp + 273.15  # Convert to Kelvin

        # Loop through all particles
        for row in self.particles:
            for particle in row:
                heat_power = 5.67037442e-8 * self.material.particle_size ** 2 * ((heat_source_temp_kelvin -
                                                                                  particle.temperature) ** 4)
                # Find center of particle, then calculate distance to center of heat source
                particle_center_x = particle.x + particle.size / 2
                particle_center_y = particle.y + particle.size / 2
                distance_px = ((heat_source_x - particle_center_x) ** 2 + (
                        heat_source_y - particle_center_y) ** 2) ** 0.5
                distance_m = distance_px * self.nm_per_pixel  # Convert to meters

                # Check if the particle is within the heat source radius
                if distance_m <= radius_nm:
                    intercepted_area = self.material.particle_size ** 2  # Convert particle size to square meters

                    # Inverse square law: Power received by the particle
                    particle_power = (heat_power * intercepted_area) / (4 * math.pi * (distance_m ** 2))

                    # Mass = Volume * Density
                    volume = self.material.particle_size ** 3  # Convert particle size to cubic meters
                    mass = volume * self.density

                    # Calculate temperature increase based on power and specific heat capacity
                    delta_temperature = particle_power * 0.05 / (self.material.specific_heat_capacity * mass)
                    particle.temperature += delta_temperature  # Add delta_temperature to the current temperature

                particle.update_color()  # Update the particle's color based on its new temperature


    # Conduction between particles
    def apply_fouriers_law(self, particle1, particle2):
        # Temperature difference
        delta_temperature = particle1.temperature - particle2.temperature

        # Temperature difference / distance between particles = temperature gradient
        temperature_gradient = delta_temperature / self.atomic_spacing

        # Same thing intercepted area from apply heat source
        area = self.material.particle_size ** 2

        # Heat transfer rate = thermal conductivity * area * temperature_gradient
        heat_transfer_rate = -self.material.thermal_conductivity * area * temperature_gradient
        return heat_transfer_rate

    def update_particle_temperatures(self):
        for i in range(self.grid_size):
            for j in range(self.grid_size):
                particle = self.particles[i][j]

                # Apply Fourier's law to each neighboring particle
                if j > 0:
                    left_neighbor = self.particles[i][j - 1]
                    heat_transfer_left = self.apply_fouriers_law(particle, left_neighbor)
                    particle.temperature += heat_transfer_left

                if j < self.grid_size - 1:
                    right_neighbor = self.particles[i][j + 1]
                    heat_transfer_right = self.apply_fouriers_law(particle, right_neighbor)
                    particle.temperature += heat_transfer_right

                if i > 0:
                    top_neighbor = self.particles[i - 1][j]
                    heat_transfer_top = self.apply_fouriers_law(particle, top_neighbor)
                    particle.temperature += heat_transfer_top

                if i < self.grid_size - 1:
                    bottom_neighbor = self.particles[i + 1][j]
                    heat_transfer_bottom = self.apply_fouriers_law(particle, bottom_neighbor)
                    particle.temperature += heat_transfer_bottom

                particle.update_color()

    # Thing to run the simulation
    def simulate_heat_transfer(self):
        self.is_simulating = True

        def update_canvas():
            if self.is_simulating:
                self.apply_heat_source()
                self.update_particle_temperatures()
                self.canvas.update()
                if not self.is_dragging:
                    # The 5 there is how many milliseconds once is this simulation run. So right now 200x per second.
                    # Can make it slower if you want
                    self.canvas.after(5, update_canvas)

        update_canvas()

    def load_background_image(self):
        # Load the image using PIL
        image = Image.open(self.material.image)

        # Resize the image to match the canvas size
        image = image.resize((self.canvas_size - 50, self.canvas_size - 50), )

        # Adjust the opacity (alpha channel) of the image
        opacity = 1  # Set the opacity level (0.0 - 1.0)
        image_with_opacity = image.copy()
        image_with_opacity.putalpha(int(255 * opacity))

        # Convert the image to a format that Tkinter can display
        self.background_image = ImageTk.PhotoImage(image_with_opacity)

        # Center the image on the canvas
        x_offset = (self.canvas_size - image_with_opacity.width) // 2
        y_offset = (self.canvas_size - image_with_opacity.height) // 2

        # Display the background image centered on the canvas
        self.canvas.create_image(x_offset, y_offset, anchor=tk.NW, image=self.background_image)

    # Creates the grid
    def create_particles_and_lines(self):
        for i in range(self.grid_size):
            row = []
            for j in range(self.grid_size):
                x = j * (self.material.particle_size + self.material.atomic_spacing) / self.nm_per_pixel + 25
                y = i * (self.material.particle_size + self.material.atomic_spacing) / self.nm_per_pixel + 25
                particle = Particle(self.canvas, x, y, self.material.particle_size / self.nm_per_pixel)
                row.append(particle)

            self.particles.append(row)

    # Move the heat source
    def move_heat_source(self, event):
        if self.is_dragging:
            self.canvas.coords(self.heat_source, event.x - 5, event.y - 5, event.x + 5, event.y + 5)

    def start_dragging(self, event):
        self.is_dragging = True

    def stop_dragging(self, event):
        if self.is_dragging:
            self.is_dragging = False
            self.simulate_heat_transfer()

    def canvas_click(self, event):
        if self.is_dragging:
            self.simulate_heat_transfer()

    def reset_simulation(self, new_material):
        self.is_simulating = False
        self.canvas.delete("all")  # Clear the canvas
        self.material = new_material  # Update the material
        self.particles = []
        # Load the new background image
        self.load_background_image()
        self.create_particles_and_lines()  # Recreate particles and lines
        self.initialize_particle_temperatures()  # Initialize particle temperatures
        self.heat_source = self.canvas.create_oval(195, 195, 205, 205, fill="red")  # Recreate the heat source
        self.canvas.tag_raise(self.heat_source)  # Make heat source in front of particles
        self.heat_source_temp = float(temperature_entry.get())  # Update heat source temperature

    def reset_simulation2(self, new_temperature):
        self.is_simulating = False
        self.canvas.delete("all")  # Clear the canvas
        self.particles = []
        # Load the new background image
        self.load_background_image()
        self.create_particles_and_lines()  # Recreate particles and lines
        self.initialize_particle_temperatures()  # Initialize particle temperatures
        self.heat_source = self.canvas.create_oval(195, 195, 205, 205, fill="red")  # Recreate the heat source
        self.canvas.tag_raise(self.heat_source)  # Make heat source in front of particles
        self.heat_source_temp = new_temperature  # Update heat source temperature

    def reset_simulation3(self, new_size):
        self.is_simulating = False
        self.canvas.delete("all")  # Clear the canvas
        self.grid_size = new_size
        self.particles = []
        self.load_background_image()
        self.nm_per_pixel = self.calculate_nm_per_pixel()
        self.create_particles_and_lines()  # Recreate particles and lines
        self.initialize_particle_temperatures()  # Initialize particle temperatures
        self.heat_source = self.canvas.create_oval(195, 195, 205, 205, fill="red")  # Recreate the heat source
        self.canvas.tag_raise(self.heat_source)  # Make heat source in front of particles
        self.heat_source_temp = float(temperature_entry.get())  # Update heat source temperature


def on_material_change(*args):
    selected_material_name = selected_material.get()
    for material in materials:
        if material.name == selected_material_name:
            simulation.reset_simulation(material)
            break


def on_unit_change(*args):
    selected_unit = selected_temperature_unit.get()
    temperature_value = float(temperature_entry.get())
    if selected_unit == "Celsius":
        new_temp = temperature_value
    elif selected_unit == "Kelvin":
        new_temp = temperature_value - 273.15
    elif selected_unit == "Fahrenheit":
        new_temp = (temperature_value - 32) * 5 / 9
    else:
        new_temp = 60

    simulation.reset_simulation2(new_temp)


def on_size_change(*args):
    new_size = int(size_entry.get())
    simulation.reset_simulation3(new_size)


root = tk.Tk()
root.title("Heat Transfer Simulation")

control_frame = tk.Frame(root)
control_frame.pack()
materials = [Material("Copper", 401.0, 385.0, 2.56e-10, 5e-8, 8960, "Copper.jpeg"),
             Material("Aluminum", 237.0, 897.0, 2.86e-10, 5e-8, 2700, "Aluminum.jpeg"),
             Material("Steel", 45.0, 500.0, 2.5e-10, 5e-8, 7850, "Steel.jpeg"),
             Material("Glass", 1.0, 840.0, 3.5e-10, 5e-8, 2500, "Glass.jpeg"),
             Material("Diamond", 2200.0, 519.0, 1.54e-10, 5e-8, 3500, "Diamond.jpeg"),
             Material("Silicon", 148.0, 710.0, 5.43e-10, 5e-8, 2330, "Silicon.jpeg"),
             Material("Platinum", 70.0, 133.0, 2.77e-10, 5e-8, 21450, "Platinum.jpeg")]
copper = materials[0]
selected_material = tk.StringVar(root)
selected_material.set(copper.name)
material_dropdown = tk.OptionMenu(control_frame, selected_material, *[material.name for material in materials])
material_dropdown.pack(side="left", padx=10)
temperature_label = tk.Label(control_frame, text="Temperature:")
temperature_label.pack(side="left", padx=0)

temperature_units = ["Celsius", "Kelvin", "Fahrenheit"]
selected_temperature_unit = tk.StringVar(root)
selected_temperature_unit.set(temperature_units[0])
temperature_unit_dropdown = tk.OptionMenu(control_frame, selected_temperature_unit,
                                          *[temperature_unit for temperature_unit in temperature_units])
temperature_unit_dropdown.pack(side="left", padx=10)


# Validation function to allow only numbers
def validate_temperature(new_value):
    if new_value == "":
        return True
    try:
        float(new_value)
        return True
    except ValueError:
        return False


# Validation function to allow only integers
def validate_size(new_value):
    if new_value == "":
        return True
    try:
        int(new_value)
        return True
    except ValueError:
        return False


# Register the validation function with the temperature entry
validate_1 = root.register(validate_temperature)

# Create the temperature entry widget with validation
temperature_entry = tk.Entry(control_frame, width=10, validate="key", validatecommand=(validate_1, "%P"))
temperature_entry.insert(0, "60")
temperature_entry.bind("<Return>", on_unit_change)
temperature_entry.pack(side="left", padx=10)

validate_2 = root.register(validate_size)
size_label = tk.Label(control_frame, text="Size:")
size_label.pack(side="left", padx=0)
size_entry = tk.Entry(control_frame, width=5, validate="key", validatecommand=(validate_2, "%P"))
size_entry.insert(0, "20")
size_entry.bind("<Return>", on_size_change)
size_entry.pack(side="left", padx=10)

move_button = tk.Button(control_frame, text="Move Heat Source")
move_button.pack(side="left", padx=10)
canvas_size = 600
grid_size = 20

simulation = HeatTransferSimulation(root, canvas_size, grid_size, copper)
move_button.bind("<ButtonPress-1>", simulation.start_dragging)

selected_material.trace("w", on_material_change)
selected_temperature_unit.trace("w", on_unit_change)

root.mainloop()
