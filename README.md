!pip install numpy scipy pandas matplotlib
import numpy as np
from scipy.integrate import solve_bvp
import pandas as pd
import matplotlib.pyplot as plt

# Cho phép hiển thị biểu đồ trong Jupyter
%matplotlib inline

# ==========================================
# 1. THIẾT LẬP THAM SỐ
# ==========================================
R_outer = 18.415e-3
R_min = 1.0e-6
h_p = 4.5e-3
h = 6.0e-6              # h=6um cho kiểm chứng
P_a = 0.1e6
P_s = 0.6e6
mu = 1.78e-5
k = 1.44e-15
K_const = (6 * k) / (h**3 * h_p)

# ==========================================
# 2. HÀM PHƯƠNG TRÌNH VÀ BIÊN
# ==========================================
def fun(r, y):
    p, dp_dr = y[0], y[1]
    term_source = (K_const / p) * (p**2 - P_s**2)
    term_geom = (1.0 / r + dp_dr / p) * dp_dr
    d2p_dr2 = term_source - term_geom
    return np.vstack((dp_dr, d2p_dr2))

def bc(ya, yb):
    return np.array([ya[1], yb[0] - P_a])

# ==========================================
# 3. CHẠY TÍNH TOÁN (Chưa vẽ ngay)
# ==========================================
mesh_resolutions = [10, 100, 1000, 10000]
r_standard_mm = np.linspace(0, R_outer*1000, 500)
df_results = pd.DataFrame({'Radius_mm': r_standard_mm})

print("--- BẮT ĐẦU TÍNH TOÁN DỮ LIỆU ---")

for N in mesh_resolutions:
    r_init = np.linspace(R_min, R_outer, N)
    p_guess = P_s - (P_s - P_a) * ((r_init - R_min) / (R_outer - R_min))**2
    dp_guess = -2 * (P_s - P_a) * (r_init - R_min) / (R_outer - R_min)**2
    y_guess = np.vstack((p_guess, dp_guess))
    
    try:
        # Tăng max_nodes lên 1 triệu để thử cứu N=10 (nhưng khả năng cao vẫn fail)
        sol = solve_bvp(fun, bc, r_init, y_guess, tol=1e-4, max_nodes=1000000)
        
        if sol.success:
            print(f"Lưới N={N}: Tính toán THÀNH CÔNG.")
            p_sol_MPa = sol.sol(r_standard_mm / 1000)[0] / 1e6
            df_results[f'Pressure_N{N}_MPa'] = p_sol_MPa
        else:
            print(f"Lưới N={N}: THẤT BẠI (Không hội tụ).")
            
    except Exception as e:
        print(f"Lưới N={N}: Lỗi code ({e})")

# Lưu dữ liệu ra Excel
df_results.to_csv('ket_qua_so_lieu.csv', index=False)
print("\n[Đã lưu dữ liệu xong. Bắt đầu vẽ hình...]")

# ==========================================
# 4. VẼ BIỂU ĐỒ (3 Riêng + 1 Gộp)
# ==========================================

# Danh sách các lưới thành công cần vẽ (Bỏ N=10 vì lỗi)
successful_meshes = [100, 1000, 10000]

# --- A. VẼ TỪNG HÌNH RIÊNG LẺ (Sửa để hiện ra màn hình) ---
for N in successful_meshes:
    col_name = f'Pressure_N{N}_MPa'
    if col_name in df_results.columns:
        plt.figure(figsize=(8, 5))
        
        plt.plot(df_results['Radius_mm'], df_results[col_name], color='blue', linewidth=2, label=f'N={N}')
        
        plt.title(f'Phân bố áp suất với lưới N={N}\n(h=6$\mu$m, Ps=0.6 MPa)')
        plt.xlabel('Bán kính r (mm)')
        plt.ylabel('Áp suất (MPa)')
        plt.axhline(P_a/1e6, color='black', linestyle=':', label='Môi trường')
        plt.grid(True)
        plt.legend()
        plt.xlim(0, 18.5)
        plt.ylim(0, 0.7)
        
        filename = f'bieu_do_rieng_N{N}.png'
        plt.savefig(filename, dpi=150)
        print(f"-> Đã lưu: {filename}")
        
        plt.show() # <--- THÊM DÒNG NÀY ĐỂ HIỆN HÌNH RA
        # plt.close() <--- XÓA HOẶC ẨN DÒNG NÀY ĐI

# --- B. VẼ HÌNH GỘP (Tất cả trong 1) ---
plt.figure(figsize=(10, 6))

for N in successful_meshes:
    col_name = f'Pressure_N{N}_MPa'
    if col_name in df_results.columns:
        # Vẽ chồng lên nhau
        lw = 2 if N == 10000 else 1.5
        ls = '-' if N == 10000 else '--'
        plt.plot(df_results['Radius_mm'], df_results[col_name], label=f'Lưới N={N}', linewidth=lw, linestyle=ls)

plt.title(f'So sánh sự hội tụ áp suất (Tổng hợp)\n(h=6$\mu$m, Ps=0.6 MPa)')
plt.xlabel('Bán kính r (mm)')
plt.ylabel('Áp suất (MPa)')
plt.axhline(P_a/1e6, color='black', linestyle=':', label='Môi trường')
plt.grid(True)
plt.legend()
plt.xlim(0, 18.5)
plt.ylim(0, 0.7)

# Lưu file gộp
plt.savefig('bieu_do_tong_hop.png', dpi=300)
print(f"-> Đã lưu: bieu_do_tong_hop.png")
plt.show() # Hiển thị hình gộp cuối cùng
