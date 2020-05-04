def dispatch_callbacks(callback_dict, callback_name, *args, **kwargs):
    callback_list = callback_dict.get(callback_name, [])
    for callback in callback_list:
        callback(*args, **kwargs)
