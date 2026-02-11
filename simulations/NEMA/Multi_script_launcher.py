from subprocess import Popen


def multi_script_launcher():
    commands = ['cd NEMA_H0 && python main.py', 'cd NEMA_H3 && python main.py', 'cd NEMA_H15 && python main.py']
    processes = [Popen(cmd, shell=True) for cmd in commands]
    for p in processes:
        p.wait()


if __name__ == '__main__':
    multi_script_launcher()
