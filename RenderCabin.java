


import java.awt.event.ComponentAdapter;
import java.awt.event.ComponentEvent;
import java.awt.event.KeyAdapter;
import java.awt.event.KeyEvent;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.event.MouseMotionListener;
import java.awt.BasicStroke;
import java.awt.Button;
import java.io.BufferedReader;
import java.io.FileReader;
import java.util.Iterator;
import java.util.Set;
import java.util.TreeSet;
import java.util.List;
import java.util.ArrayList;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.Image;
import java.awt.Panel;
import java.awt.Polygon;
import java.awt.TextField;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.Canvas;
import java.awt.Checkbox;
import java.awt.Color;
import java.awt.Frame;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.awt.image.BufferedImage;
import java.awt.image.MemoryImageSource;


// Main class of the project
@SuppressWarnings("serial")
public class RenderCabin extends Frame implements ActionListener {
TextField text_f1, text_f2, text_f3;
	RenderCanvas renderer;
	Checkbox choice;


	// Constructor
	public RenderCabin(String file, int yaw, int pitch, double focal, double size) {
		super("Cabin Model");
		try {
			Luminance.yaw = yaw * Math.PI / 180;
			Luminance.pitch = pitch * Math.PI / 180;
			Luminance.focal = focal;
			Luminance.size = size;
			BufferedReader br = new BufferedReader(new FileReader(file));
			int num = Integer.parseInt(br.readLine());
			Vertex vertices[] = new Vertex[num];
			for (int i = 0; i < num; i++) {
				vertices[i] = Vertex.parseVertex(br.readLine());
				vertices[i].projectToScreen();
			}
			Polygon3D.vertices = vertices;
			num = Integer.parseInt(br.readLine());
			Polygon3D polygons[] = new Polygon3D[num];
			for (int i = 0; i < num; i++)
				polygons[i] = Polygon3D.parsePolygon3D(br.readLine());
			renderer = new RenderCanvas(vertices, polygons);
			br.close();
		} catch (Exception e) {
			System.out.println(e);
		}
		Panel controls = new Panel();
		choice = new Checkbox("Convex", null, false);
		controls.add(choice);
		
		Button button = new Button("Wireframe");
		button.addActionListener(this);
		controls.add(button);
		button = new Button("Polygon Color");
		button.addActionListener(this);
		controls.add(button);
		button = new Button("Constant Shading");
		button.addActionListener(this);
		controls.add(button);
		button = new Button("Perpixel Shading");
		button.addActionListener(this);
		controls.add(button);
		// add canvas and control panel
		add("Center", renderer);
		add("South", controls);
		addWindowListener(new ExitListener());
		//creating text field for lightening with just one column
				text_f1 = new TextField(Double.toString(Luminance.c1), 1);
			controls.add(text_f1);
			text_f2 = new TextField(Double.toString(Luminance.c2), 1);
				controls.add(text_f2);
				text_f3 = new TextField(Double.toString(Luminance.c3), 1);
				controls.add(text_f3);
	}

	class ExitListener extends WindowAdapter {
		public void windowClosing(WindowEvent e) {
			System.exit(0);
		}
	}
	
	// Action listener for buttons
	public void actionPerformed(ActionEvent e) {
		renderer.withHole = !choice.getState();
		
		if (((Button) e.getSource()).getLabel().equals("Wireframe"))
			renderer.rendering(0);
		else if (((Button) e.getSource()).getLabel().equals("Polygon Color"))
			renderer.rendering(1);
		else if (((Button) e.getSource()).getLabel().equals("Constant Shading"))
			renderer.rendering(2);
		else if (((Button) e.getSource()).getLabel().equals("Perpixel Shading"))
			renderer.rendering(3);
		renderer.repaint();
Luminance.c1 = Double.parseDouble(text_f1.getText());
Luminance.c2 = Double.parseDouble(text_f2.getText());
Luminance.c3 = Double.parseDouble(text_f3.getText());
	}

	public static void main(String[] args) {
		int raw = 200, pitch = 20;
		if (args.length >= 1)
			raw = Integer.parseInt(args[0]);
		if (args.length >= 2)
			pitch = Integer.parseInt(args[1]);
		RenderCabin window = new RenderCabin("cabin.mdl", raw, pitch, 15, 20);
		window.setSize(600, 600);
		window.setVisible(true);
	}
}
@SuppressWarnings("serial")
class RenderCanvas extends Canvas {
	// image for display.
	int commonrendering;
	Image image;
	boolean withHole = true;
	
	// the triangles for rendering.
	Vertex vertices[];
	Polygon3D polygons[];
	// canvas size
	int width, height;

	// initialize the canvas
	public RenderCanvas(Vertex ver[], Polygon3D poly[]) {
		Polygon3D.vertices = vertices = ver;
		polygons = poly;
		addComponentListener(new ResizeListener());
		commonrendering = 0;
		DragListener drag = new DragListener();
		addMouseListener(drag);
		addMouseMotionListener(drag);
		addKeyListener(new ArrowListener());
		addComponentListener(new ResizeListener());
	}

	class ResizeListener extends ComponentAdapter {
		public void componentResized(ComponentEvent e) {
			width = getWidth();
			height = getHeight();
			rendering(-1);
		}
	}

	// Action listener for mouse
	class DragListener extends MouseAdapter implements MouseMotionListener {
		int ending_X, ending_Y;

		public void mousePressed(MouseEvent e) {
			ending_X = e.getX();
			ending_Y = e.getY();
		}

		// update pitch and yaw angles when the mouse is dragged.
		public void mouseDragged(MouseEvent e) {
			Luminance.yaw -= (e.getX() - ending_X) * Math.PI / 180;
			Luminance.pitch -= (e.getY() - ending_Y) * Math.PI / 180;
			ending_X = e.getX();
			ending_Y = e.getY();
			reproject();
			rendering(-1);
			repaint();
		}
	}
	public void reproject() {
		for (int j = 0; j < Polygon3D.vertices.length; j++) {
			Polygon3D.vertices[j].projectToScreen();
		}
	}

	// Action listener for keyboard
	class ArrowListener extends KeyAdapter {
		public void keyPressed(KeyEvent e) {
			if (e.getKeyCode() == KeyEvent.VK_DOWN && Luminance.focal > 5)
				Luminance.focal--;
			else if (e.getKeyCode() == KeyEvent.VK_UP && Luminance.focal < 30)
				Luminance.focal++;
			else if (e.getKeyCode() == KeyEvent.VK_LEFT && Luminance.size > 10)
				Luminance.size--;
			else if (e.getKeyCode() == KeyEvent.VK_RIGHT && Luminance.size < 30)
				Luminance.size++;
			reproject();
			rendering(-1);
			repaint();
		}
	}

	
	public void rendering(int choice) {
		if (choice == 0) {
			commonrendering = choice;
			renderWireframe();
		} else if (choice == 1) {
			commonrendering = choice;
			Rendering(false, false);
		} else if (choice == 2) {
			commonrendering = choice;
			Rendering(true, false);
		} else if (choice == 3) {
			commonrendering = choice;
			Rendering(false, true);
		} else {
			rendering(commonrendering);
		}
	}

	// render the wireframe of triangle mesh
	public void renderWireframe() {
		BufferedImage srcImage = new BufferedImage(width, height, BufferedImage.TYPE_INT_RGB);
		Graphics2D g2D = srcImage.createGraphics();
		g2D.setBackground(Color.cyan);
		g2D.clearRect(0, 0, width, height);
		g2D.setStroke(new BasicStroke(4));
		for (int i = 0; i < polygons.length; i++) {
			g2D.setColor(polygons[i].color);
			g2D.draw(polygons[i].toPolygon2D(width, height, withHole));
		}
		image = srcImage;
	} 
// determine the polygon intersections
	public Get_intersection polygonintersections(int current_y, Polygon3D polygon3D, boolean diff_color) {
		Polygon polygon = polygon3D.toPolygon2D(width, height, withHole);
		polygon3D.Normal_ConstantColor(width, height);
		Get_intersection intersections = new Get_intersection();
		List<Edge> edges = new ArrayList<>();
		List<IgnoreV> ignoreV = new ArrayList<>();
		for (int i = 0; i < polygon.npoints; i++) {
			int next = (i + 1) % polygon.npoints;
			int y_p1 = polygon.ypoints[i];
			int y_p2 = polygon.ypoints[next];
			int previous = ((i - 1) + polygon.npoints) % polygon.npoints;
			int y_p0 = polygon.ypoints[previous];

			Vector3D v_1 = new Vector3D(polygon.xpoints[i], polygon.ypoints[i], Polygon3D.vertices[polygon3D.id[i]].z);
			Vector3D v_2 = new Vector3D(polygon.xpoints[next], polygon.ypoints[next], Polygon3D.vertices[polygon3D.id[next]].z);
			Edge _edge = new Edge(v_1, v_2, i, next);

			
			if (y_p1 != y_p2 && current_y <= _edge.start.y && current_y > _edge.end.y) {
				edges.add(_edge);
			}

			if ((y_p0 - y_p1) * (y_p2 - y_p1) > 0) {
				ignoreV.add(new IgnoreV(v_1, i));
			}
		}
		Iterator<Edge> eiter = edges.iterator();
		while (eiter.hasNext()) {
			Edge _edge = eiter.next();
			Vector3D _current = _edge.current(current_y);
			int indexing = _edge.indexing(current_y);
			boolean ignore = false;
            for (int j = 0; j < ignoreV.size(); j++)
            {
                if (ignoreV.get(j).compareTo(new IgnoreV(_current,indexing))==0)
                {
                    ignore = true;
                }
            }
            if(current_y != _edge.start.y || !ignore)
			if (diff_color) {
				intersections.VertexAdding((int) _current.x, _current.z, polygon3D, polygon3D.ncolor, indexing);
			} else {
				intersections.VertexAdding((int) _current.x, _current.z, polygon3D, polygon3D.color, indexing);
			}
		}
		return intersections;
	}

	public void Rendering(boolean changeColor, boolean perpixel) {
		int _pixels[] = new int[width * height];
		for (int j = 0; j < height; j++) {
			Get_intersection[] inter = new Get_intersection[polygons.length];
			for (int l = 0; l < polygons.length; l++) {
				inter[l] = polygonintersections(j, polygons[l], changeColor);
			}
			for (int i = 0; i < width; i++) {
				Color c = Color.cyan;
				Polygon3D poly3d = null;
				double zmin = Double.POSITIVE_INFINITY;
				for (int p = 0; p < polygons.length; p++) {
					double z = inter[p].CheckInside(i);
					if (!Double.isNaN(z) && z < zmin) {
						zmin = z;
						c = inter[p].color;
						poly3d = inter[p].polygon3d;
					}
				}
				if (perpixel && poly3d != null) {
					Vector3D wpoint = new Vector3D(i, j, zmin).backToWorld(width, height);
					_pixels[j * width + i] = Luminance.PhongModel(wpoint, poly3d.normal, c, width, height).getRGB();
				} else {
					_pixels[j * width + i] = c.getRGB();
				}
			}
		}
		image = createImage(new MemoryImageSource(width, height, _pixels, 0, width));
	}

	public void paint(Graphics g) {
		// draw the rendering result.
		g.drawImage(image, 0, 0, this);
		// display current parameter values
		g.setColor(Color.black);
		g.drawString("Pitch = " + (int) (Luminance.pitch * 180 / Math.PI), 10, 20);
		g.drawString("Yaw = " + (int) (Luminance.yaw * 180 / Math.PI), 10, 40);
		g.drawString("Focal = " + Luminance.focal, 10, 60);
		g.drawString("Size = " + Luminance.size, 10, 80);
	}
}

// Vertex definition: position and texture coords
class Vertex {
	double start_x, start_y, start_z;
	double x, y, z;
	double u, v;

	public Vertex(double V_x, double V_y, double V_z) {
		start_x = V_x;
		start_y = V_y;
		start_z = V_z;
	}

	// parse a vertex from a string
	public static Vertex parseVertex(String input) {
		String tokens[] = input.split("\\s+");
		Vertex vert = new Vertex(Double.parseDouble(tokens[0]), Double.parseDouble(tokens[1]),
				Double.parseDouble(tokens[2]));
		vert.u = Double.parseDouble(tokens[3]);
		vert.v = Double.parseDouble(tokens[4]);
		return vert;
	}

	// project the vertex to normalized device coordinate
	public void projectToScreen() {
		// rotate about y-axis
		double x1 = Math.cos(Luminance.yaw) * start_x - Math.sin(Luminance.yaw) * start_z;
		double z1 = Math.sin(Luminance.yaw) * start_x + Math.cos(Luminance.yaw) * start_z;
		// rotate about x-axis
		double y2 = Math.cos(Luminance.pitch) * start_y + Math.sin(Luminance.pitch) * z1;
		double z2 = -Math.sin(Luminance.pitch) * start_y + Math.cos(Luminance.pitch) * z1;
		// translate object to focal plane
		z2 += Luminance.focal;
		// perspective projection
		x = x1 * 2 * Luminance.focal / Luminance.size / z2;
		y = y2 * 2 * Luminance.focal / Luminance.size / z2;
		z = 1 - 2 / z2;
	}

}

// Polygon definition
class Polygon3D {
	static Vertex vertices[];

	// indicies of all vertices
	int id[];
	Color color;
	Color ncolor;

	Vector3D normal;

	public Polygon3D(int num) {
		id = new int[num];
	}

	// parse a 3D polygon from a string
	public static Polygon3D parsePolygon3D(String input) {
		String tokens[] = input.split("\\s");
		int num = Integer.parseInt(tokens[0]);
		Polygon3D poly = new Polygon3D(num);
		poly.color = new Color(Float.parseFloat(tokens[1]), Float.parseFloat(tokens[2]), Float.parseFloat(tokens[3]));
		for (int n = 0; n < num; n++)
			poly.id[n] = Integer.parseInt(tokens[n + 4]);
		return poly;
	}

	// convert 3D polygon to a 2D polygon
	public Polygon toPolygon2D(int width, int height, boolean hole) {
		Polygon poly = new Polygon();
		double scale = Math.min(width, height);
		int num = hole ? id.length : Math.min(4, id.length);
		for (int n = 0; n < num; n++)
			poly.addPoint((int) (vertices[id[n]].x * scale + width / 2),
					(int) (-vertices[id[n]].y * scale + height * 3 / 4));
		return poly;
	}
// it will compute Normal vector and Constant color for our polygon
	public void Normal_ConstantColor(int width, int height) {

		Vertex v2D_0 = vertices[id[0]];
		Vertex v2D_1 = vertices[id[1]];
		Vertex v2D_2 = vertices[id[2]];

		Vector3D v3D_01 = new Vector3D(v2D_1.start_x - v2D_0.start_x, v2D_1.start_y - v2D_0.start_y, v2D_1.start_z - v2D_0.start_z);
		Vector3D v3D_02 = new Vector3D(v2D_2.start_x - v2D_0.start_x, v2D_2.start_y - v2D_0.start_y, v2D_2.start_z - v2D_0.start_z);

		Vector3D n = new Vector3D(v3D_01.y * v3D_02.z - v3D_01.z * v3D_02.y, v3D_01.z * v3D_02.x - v3D_01.x * v3D_02.z,
				v3D_01.x * v3D_02.y - v3D_01.y * v3D_02.x);

		Vector3D L = Luminance.view(new Vector3D(v2D_0.start_x, v2D_0.start_y, v2D_0.start_z), width, height);

		if (L.x * n.x + L.y * n.y + L.z * n.z >= 0) {
			normal = n;
		} else {
			normal = new Vector3D(-n.x, -n.y, -n.z);
		}

		double rx = 0, ry = 0, rz = 0;
		for (int rr = 0; rr < id.length; rr++) {
			rx += vertices[id[rr]].x;
			ry += vertices[id[rr]].y;
			rz += vertices[id[rr]].z;
		}
		rx = rx / id.length;
		ry = ry / id.length;
		rz = rz / id.length;

		ncolor = Luminance.PhongModel(new Vector3D(rx, ry, rz), normal, color, width, height);
	}
}

class Edge {
	Vector3D start;
	Vector3D end;
	int start_i;
	int end_i;

	public Edge(Vector3D _v1, Vector3D v_2, int _i1, int _i2) {
		if (_v1.y > v_2.y) {
			start = _v1;
			end = v_2;
			start_i = _i1;
			end_i = _i2;
		} else {
			start = v_2;
			end = _v1;
			start_i = _i2;
			end_i = _i1;
		}
	}

	public Vector3D current(int cur_y) {
		double k = (double) (cur_y - end.y) / (double) (start.y - end.y);
		double x = k * (start.x - end.x) + end.x;
		double z = k * (start.z - end.z) + end.z;
		return new Vector3D(x, cur_y, z);
	}

	public int indexing(int current_y) {
		double k = (start.y + end.y) / 2;
		if (current_y >= k) {
			return start_i;
		} else {
			return end_i;
		}
	}
}

class IgnoreV implements Comparable<IgnoreV> {
	Vector3D v;
	int index;

	public IgnoreV(Vector3D _v, int _index) {
		v = _v;
		index = _index;
	}

	public int compareTo(IgnoreV o) {
		double compareto_x = v.x - o.v.x;
		double compareto_y = v.y - o.v.y;
		double compareto_z = v.z - o.v.z;
		int c_i = index - o.index;
		return (int) (compareto_x * compareto_x + compareto_y * compareto_y + compareto_z * compareto_z + c_i * c_i);
	}
}

class Get_intersection {

	Set<Common> points;

	public class Common implements Comparable<Common> {
		int index;double z;int x;

		public Common(int _x, double _z, int _index) {
			x = _x;
			z = _z;
			index = _index;
		}

		public int compareTo(Common o) {
			int c1 = x - o.x;
			return (c1 == 0) ? (index - o.index) : c1;
		}

		@Override
		public int hashCode() {
			return (31 * index + x);
		}

		@Override
		public boolean equals(Object o) {
			if (!(o instanceof Common))
				return false;
			Common object = (Common) o;
			if (x == object.x) {
				return index == object.index;
			}
			return false;
		}
	}
	Color color;
	Polygon3D polygon3d;


	public Get_intersection() {
		color = Color.cyan;
		points = new TreeSet<>();
	}

	public void VertexAdding(int x, double z, Polygon3D p, Color c, int _index) {
		polygon3d = p;
		if (!color.equals(c)) {
			color = c;
		}
		if (!points.contains(new Common(x, z, _index))) {
			points.add(new Common(x, z, _index));
		} else {
			points.remove(new Common(x, z, _index));
		}
	}
// it cheack if the point is inside or not
	public double CheckInside(int x) {
		Iterator<Common> sit = points.iterator();
		while (sit.hasNext()) {
			Common n1 = sit.next();
			if (sit.hasNext()) {
				Common n2 = sit.next();
				if (x >= n1.x && x <= n2.x) {
					if (x == n1.x && x == n2.x) {
						return n1.z;
					} else {
						double z1 = n1.z;
						double z2 = n2.z;
						double k = (double) (x - n1.x) / (double) (n2.x - n1.x);
						return k * (z2 - z1) + z1;
					}
				}
			}
		}
		return Double.NaN;
	}
}

class Vector3D {
	double x, y, z;

	// constructors
	public Vector3D(double a, double b, double c) {
		x = a;
		y = b;
		z = c;
	}
// make normalization
	public void getnorm() {
		double norm = Math.sqrt(x * x + y * y + z * z);
		x = x / norm;
		y = y / norm;
		z = z / norm;
	}

	public Vector3D backToWorld(int width, int height) {

		double scale = Math.min(width, height);
		double sx = (x - width / 2) / scale;
		double sy = (y - height * 3 / 4) / scale;

		double z2 = 2 / (double) (1 - z);
		double y2 = sy * z2 * Luminance.size / Luminance.focal / 2;
		double x2 = sx * z2 * Luminance.size / Luminance.focal / 2;

		z2 -= Luminance.focal;

		double y_1 = Math.cos(Luminance.pitch) * y2 - Math.sin(Luminance.pitch) * z2;
		double z_1 = Math.sin(Luminance.pitch) * y2 + Math.cos(Luminance.pitch) * z2;

		double x_0 = Math.cos(Luminance.yaw) * x2 + Math.sin(Luminance.yaw) * z_1;
		double z_0 = -Math.sin(Luminance.yaw) * x2 + Math.cos(Luminance.yaw) * z_1;

		return new Vector3D(x_0, y_1, z_0);
	}
}

class Luminance {
	static int n = 5;
	static Vector3D source = new Vector3D(-2, 3, 4);
	static Vector3D intensity = new Vector3D(1, 1, 1);
	static double yaw, pitch, focal, size;
	static double c1 = 1;
	static double c2 = 0;
	static double c3 = 0.01;
	static double kd = 0.7;
	static double ks = 0.3;



	static public Vector3D intensity(Vector3D point) {
		double d = Math.sqrt(Math.pow((point.x - source.x), 2) + Math.pow((point.y - source.y), 2)
				+ Math.pow((point.z - source.z), 2));
		double dis = Math.min(1, 1 / (c1 + c2 * d + c3 * d * d));
		return new Vector3D(dis * intensity.x, dis * intensity.y, dis * intensity.z);
	}

	static public Vector3D reflectVector(Vector3D a, Vector3D b) {
		double norm = Math.sqrt(b.x * b.x + b.y * b.y + b.z * b.z);
		Vector3D u = new Vector3D(b.x / norm, b.y / norm, b.z / norm);
		double t = a.x * u.x + a.y * u.y + a.z * u.z;
		Vector3D d = new Vector3D(t * u.x, t * u.y, t * u.z);
		return new Vector3D(2 * d.x - a.x, 2 * d.y - a.y, 2 * d.z - a.z);
	}

	static public Vector3D view(Vector3D P, int width, int height) {
		Vector3D v = new Vector3D(0, 0, 0).backToWorld(width, height);
		return new Vector3D(v.x - P.x, v.y - P.y, v.z - P.z);
	}

	static public Color PhongModel(Vector3D P, Vector3D N, Color c, int width, int height) {
		N.getnorm();
		Vector3D L = new Vector3D(source.x - P.x, source.y, source.z - P.z);
		L.getnorm();
		Vector3D R = reflectVector(L, N);
		R.getnorm();
		Vector3D I = intensity(P);
		Vector3D V = view(P, width, height);
		V.getnorm();
		double diffuse = kd * (N.x * L.x + N.y * L.y + N.z * L.z);
		Color dcolor = new Color(pin(I.x * c.getRed() * diffuse), pin(I.y * c.getGreen() * diffuse),
				pin(I.z * c.getBlue() * diffuse));
		double specular = ks * (V.x * R.x + V.y * R.y + V.z * R.z);
		Color scolor = new Color(pin(I.x * c.getRed() * specular), pin(I.y * c.getGreen() * specular),
				pin(I.z * c.getBlue() * specular));
		return new Color(pin(dcolor.getRed() + scolor.getRed()), pin(dcolor.getGreen() + scolor.getGreen()),
				pin(dcolor.getBlue() + scolor.getBlue()));
	}

	static public int pin(double value) {
		int v = (int) value;
		v = Math.min(v, 255);
		return Math.max(v, 0);
	}
}