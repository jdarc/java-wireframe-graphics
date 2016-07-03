import java.awt.*;
import java.awt.geom.Line2D;

public class Graphics3D {

    public enum Mode {
        LINES,
        POLYGON
    }

    private final Graphics2D graphics;
    private final float[] world = new float[16];
    private final float[] view = new float[16];
    private final float[] proj = new float[16];
    private final float[] comb = new float[16];
    private final Vertex[] vertices = new Vertex[256];
    private final Rectangle viewport = new Rectangle(0, 0, 0, 0);
    private Mode drawMode = Mode.LINES;
    private boolean dirty = true;
    private int count = 0;

    public Graphics3D(Graphics2D graphics) {
        this.graphics = graphics;
        for (int i = 0; i < vertices.length; i++) {
            vertices[i] = new Vertex();
        }
    }

    public Color getColor() {
        return graphics.getColor();
    }

    public void setColor(Color color) {
        graphics.setColor(color);
    }

    public void viewport(int x, int y, int width, int height) {
        viewport.setBounds(x, y, width, height);
    }

    public void clear() {
        graphics.fillRect(viewport.x, viewport.y, viewport.width, viewport.height);
    }

    public void identity() {
        dirty = true;
        makeIdentity(world);
    }

    public void scale(float x, float y, float z) {
        dirty = true;
        float[] transform = makeIdentity(new float[16]);
        transform[0x0] = x;
        transform[0x5] = y;
        transform[0xA] = z;
        concatenate(transform, world, world);
    }

    public void rotate(float angle, float x, float y, float z) {
        dirty = true;
        float len = (float)Math.sqrt(x * x + y * y + z * z);
        float nx = x / len;
        float ny = y / len;
        float nz = z / len;
        float rad = -angle;
        float c = (float)Math.cos(rad);
        float s = (float)Math.sin(rad);
        float t = 1f - c;
        float[] transform = makeIdentity(new float[16]);
        transform[0x0] = c + nx * nx * t;
        transform[0x1] = nx * ny * t - nz * s;
        transform[0x2] = nx * nz * t + ny * s;
        transform[0x4] = nx * ny * t + nz * s;
        transform[0x5] = c + ny * ny * t;
        transform[0x6] = ny * nz * t - nx * s;
        transform[0x8] = nx * nz * t - ny * s;
        transform[0x9] = ny * nz * t + nx * s;
        transform[0xA] = c + nz * nz * t;
        concatenate(transform, world, world);
    }

    public void translate(float x, float y, float z) {
        dirty = true;
        float[] transform = makeIdentity(new float[16]);
        transform[0xC] = x;
        transform[0xD] = y;
        transform[0xE] = z;
        concatenate(transform, world, world);
    }

    public void view(float eyex, float eyey, float eyez, float centerx, float centery, float centerz, float upx, float upy, float upz) {
        dirty = true;
        float z0 = eyex - centerx;
        float z1 = eyey - centery;
        float z2 = eyez - centerz;
        float lz = (float)(1.0 / Math.sqrt(z0 * z0 + z1 * z1 + z2 * z2));
        float x0 = upy * z2 - upz * z1;
        float x1 = upz * z0 - upx * z2;
        float x2 = upx * z1 - upy * z0;
        float lx = (float)(1.0 / Math.sqrt(x0 * x0 + x1 * x1 + x2 * x2));
        float y0 = z1 * x2 - z2 * x1;
        float y1 = z2 * x0 - z0 * x2;
        float y2 = z0 * x1 - z1 * x0;
        float ly = (float)(1.0 / Math.sqrt(y0 * y0 + y1 * y1 + y2 * y2));
        view[0x0] = x0 * lx;
        view[0x1] = y0 * ly;
        view[0x2] = z0 * lz;
        view[0x4] = x1 * lx;
        view[0x5] = y1 * ly;
        view[0x6] = z1 * lz;
        view[0x8] = x2 * lx;
        view[0x9] = y2 * ly;
        view[0xA] = z2 * lz;
        view[0xC] = -(x0 * eyex + x1 * eyey + x2 * eyez) * lx;
        view[0xD] = -(y0 * eyex + y1 * eyey + y2 * eyez) * ly;
        view[0xE] = -(z0 * eyex + z1 * eyey + z2 * eyez) * lz;
        view[0xF] = 1f;
    }

    public void orthographic(float left, float right, float bottom, float top, float near, float far) {
        dirty = true;
        proj[0x0] = 2f / (right - left);
        proj[0x5] = 2f / (top - bottom);
        proj[0xA] = -2f / (far - near);
        proj[0xB] = 0f;
        proj[0xC] = -(right + left) / (right - left);
        proj[0xD] = -(top + bottom) / (top - bottom);
        proj[0xE] = -(far + near) / (far - near);
        proj[0xF] = 1f;
    }

    public void perspective(float fovy, float aspect, float near, float far) {
        dirty = true;
        double tan = Math.tan(fovy / 2.0);
        proj[0x0] = (float)(1.0 / (tan * aspect));
        proj[0x5] = (float)(1.0 / tan);
        proj[0xA] = (far + near) / (near - far);
        proj[0xB] = -1f;
        proj[0xC] = 0f;
        proj[0xD] = 0f;
        proj[0xE] = 2f * far * near / (near - far);
        proj[0xF] = 0f;
    }

    public void begin(Mode mode) {
        drawMode = mode;
        count = 0;
    }

    public void vertex(float x, float y, float z) {
        vertices[count++].copy(x, y, z, 0f);
    }

    public void end() {
        Color oldcolor = getColor();
        float[] matrix = updateTransform();
        if (drawMode == Mode.LINES) {
            if (count > 1) {
                for (int i = 0; i < count - 1; i += 2) {
                    drawLine(vertices[i].project(matrix), vertices[i + 1].project(matrix));
                }
            }
        } else if (drawMode == Mode.POLYGON) {
            if (count > 2) {
                for (int i = 0; i < count; i++) {
                    vertices[i].project(matrix);
                }
                if (isBackFacing(vertices[0], vertices[1], vertices[2])) {
                    setColor(new Color(oldcolor.getRed(), oldcolor.getGreen(), oldcolor.getBlue(), 32));
                }
                int pad = (count << 1) - 1;
                Vertex v0 = vertices[pad - 1];
                Vertex v1 = vertices[pad];
                drawLine(v0.copy(vertices[pad >> 1]), v1.copy(vertices[0]));
                for (int i = pad - 2; i > 0; i -= 2) {
                    drawLine(v0.copy(vertices[(i + 1 >> 1) - 1]), v1.copy(vertices[i + 1 >> 1]));
                }
            }
        }
        setColor(oldcolor);
    }

    private boolean isBackFacing(Vertex a, Vertex b, Vertex c) {
        return (c.x / c.w - a.x / a.w) * (a.y / a.w - b.y / b.w) < (a.y / a.w - c.y / c.w) * (b.x / b.w - a.x / a.w);
    }

    private void drawLine(Vertex a, Vertex b) {
        if (clip(a, b)) {
            a.screen(viewport);
            b.screen(viewport);
            graphics.draw(new Line2D.Float(a.x, a.y, b.x, b.y));
        }
    }

    private float[] updateTransform() {
        if (dirty) {
            dirty = false;
            concatenate(proj, view, comb);
            concatenate(comb, world, comb);
        }
        return comb;
    }

    private float[] makeIdentity(float[] dst) {
        dst[0x0] = 1f; dst[0x4] = 0f; dst[0x8] = 0f; dst[0xC] = 0f;
        dst[0x1] = 0f; dst[0x5] = 1f; dst[0x9] = 0f; dst[0xD] = 0f;
        dst[0x2] = 0f; dst[0x6] = 0f; dst[0xA] = 1f; dst[0xE] = 0f;
        dst[0x3] = 0f; dst[0x7] = 0f; dst[0xB] = 0f; dst[0xF] = 1f;
        return dst;
    }

    private float[] concatenate(float[] a, float[] b, float[] dst) {
        float m0 = a[0x0] * b[0x0] + a[0x4] * b[0x1] + a[0x8] * b[0x2] + a[0xC] * b[0x3];
        float m1 = a[0x1] * b[0x0] + a[0x5] * b[0x1] + a[0x9] * b[0x2] + a[0xD] * b[0x3];
        float m2 = a[0x2] * b[0x0] + a[0x6] * b[0x1] + a[0xA] * b[0x2] + a[0xE] * b[0x3];
        float m3 = a[0x3] * b[0x0] + a[0x7] * b[0x1] + a[0xB] * b[0x2] + a[0xF] * b[0x3];
        float m4 = a[0x0] * b[0x4] + a[0x4] * b[0x5] + a[0x8] * b[0x6] + a[0xC] * b[0x7];
        float m5 = a[0x1] * b[0x4] + a[0x5] * b[0x5] + a[0x9] * b[0x6] + a[0xD] * b[0x7];
        float m6 = a[0x2] * b[0x4] + a[0x6] * b[0x5] + a[0xA] * b[0x6] + a[0xE] * b[0x7];
        float m7 = a[0x3] * b[0x4] + a[0x7] * b[0x5] + a[0xB] * b[0x6] + a[0xF] * b[0x7];
        float m8 = a[0x0] * b[0x8] + a[0x4] * b[0x9] + a[0x8] * b[0xA] + a[0xC] * b[0xB];
        float m9 = a[0x1] * b[0x8] + a[0x5] * b[0x9] + a[0x9] * b[0xA] + a[0xD] * b[0xB];
        float ma = a[0x2] * b[0x8] + a[0x6] * b[0x9] + a[0xA] * b[0xA] + a[0xE] * b[0xB];
        float mb = a[0x3] * b[0x8] + a[0x7] * b[0x9] + a[0xB] * b[0xA] + a[0xF] * b[0xB];
        float mc = a[0x0] * b[0xC] + a[0x4] * b[0xD] + a[0x8] * b[0xE] + a[0xC] * b[0xF];
        float md = a[0x1] * b[0xC] + a[0x5] * b[0xD] + a[0x9] * b[0xE] + a[0xD] * b[0xF];
        float me = a[0x2] * b[0xC] + a[0x6] * b[0xD] + a[0xA] * b[0xE] + a[0xE] * b[0xF];
        float mf = a[0x3] * b[0xC] + a[0x7] * b[0xD] + a[0xB] * b[0xE] + a[0xF] * b[0xF];
        dst[0x0] = m0; dst[0x1] = m1; dst[0x2] = m2; dst[0x3] = m3;
        dst[0x4] = m4; dst[0x5] = m5; dst[0x6] = m6; dst[0x7] = m7;
        dst[0x8] = m8; dst[0x9] = m9; dst[0xA] = ma; dst[0xB] = mb;
        dst[0xC] = mc; dst[0xD] = md; dst[0xE] = me; dst[0xF] = mf;
        return dst;
    }

    private boolean clip(Vertex a, Vertex b) {
        for (int i = 0; i < 8; i++) {
            int amask = (a.x > +a.w ? 0x20 : 0x00) | (a.y > +a.w ? 0x10 : 0x00) | (a.z > +a.w ? 0x08 : 0x00) |
                        (a.x < -a.w ? 0x04 : 0x00) | (a.y < -a.w ? 0x02 : 0x00) | (a.z < -a.w ? 0x01 : 0x00);
            int bmask = (b.x > +b.w ? 0x20 : 0x00) | (b.y > +b.w ? 0x10 : 0x00) | (b.z > +b.w ? 0x08 : 0x00) |
                        (b.x < -b.w ? 0x04 : 0x00) | (b.y < -b.w ? 0x02 : 0x00) | (b.z < -b.w ? 0x01 : 0x00);
            if (amask == 0 && bmask == 0) return true;
            else if ((amask & bmask) != 0) return false;
            else if ((amask & 0x20) > 0) a.lerp(a, b, (a.w - a.x) / (a.w - a.x - b.w + b.x));
            else if ((amask & 0x10) > 0) a.lerp(a, b, (a.w - a.y) / (a.w - a.y - b.w + b.y));
            else if ((amask & 0x08) > 0) a.lerp(a, b, (a.w - a.z) / (a.w - a.z - b.w + b.z));
            else if ((amask & 0x04) > 0) a.lerp(a, b, (a.w + a.x) / (a.w + a.x - b.w - b.x));
            else if ((amask & 0x02) > 0) a.lerp(a, b, (a.w + a.y) / (a.w + a.y - b.w - b.y));
            else if ((amask & 0x01) > 0) a.lerp(a, b, (a.w + a.z) / (a.w + a.z - b.w - b.z));
            else if ((bmask & 0x20) > 0) b.lerp(b, a, (b.w - b.x) / (b.w - b.x - a.w + a.x));
            else if ((bmask & 0x10) > 0) b.lerp(b, a, (b.w - b.y) / (b.w - b.y - a.w + a.y));
            else if ((bmask & 0x08) > 0) b.lerp(b, a, (b.w - b.z) / (b.w - b.z - a.w + a.z));
            else if ((bmask & 0x04) > 0) b.lerp(b, a, (b.w + b.x) / (b.w + b.x - a.w - a.x));
            else if ((bmask & 0x02) > 0) b.lerp(b, a, (b.w + b.y) / (b.w + b.y - a.w - a.y));
            else if ((bmask & 0x01) > 0) b.lerp(b, a, (b.w + b.z) / (b.w + b.z - a.w - a.z));
        }
        return false;
    }

    private class Vertex {

        float x = 0f;
        float y = 0f;
        float z = 0f;
        float w = 0f;

        Vertex copy(float x, float y, float z, float w) {
            this.x = x;
            this.y = y;
            this.z = z;
            this.w = w;
            return this;
        }

        Vertex copy(Vertex other) {
            return copy(other.x, other.y, other.z, other.w);
        }

        Vertex screen(Rectangle viewport) {
            return copy((1f + x / w) * 0.5f * viewport.width + viewport.x,
                        (1f - y / w) * 0.5f * viewport.height + viewport.y,
                        (1f + z / w) * 0.5f, 1f);
        }

        Vertex project(float[] mat) {
            return copy(mat[0x0] * x + mat[0x4] * y + mat[0x8] * z + mat[0xC],
                        mat[0x1] * x + mat[0x5] * y + mat[0x9] * z + mat[0xD],
                        mat[0x2] * x + mat[0x6] * y + mat[0xA] * z + mat[0xE],
                        mat[0x3] * x + mat[0x7] * y + mat[0xB] * z + mat[0xF]);
        }

        Vertex lerp(Vertex from, Vertex to, Float t) {
            return copy(from.x + t * (to.x - from.x),
                        from.y + t * (to.y - from.y),
                        from.z + t * (to.z - from.z),
                        from.w + t * (to.w - from.w));
        }
    }
}
