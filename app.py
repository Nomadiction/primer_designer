import os
import threading
import webview
import time

def start_streamlit():
    os.system("streamlit run frontend/gui_app.py --server.headless true --server.port 8501")

threading.Thread(target=start_streamlit, daemon=True).start()

time.sleep(2.5)  

webview.create_window("Cloning Assistant", "http://localhost:8501", width=1200, height=800)
webview.start()
