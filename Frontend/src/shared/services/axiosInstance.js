import axios from 'axios';
import { BACKEND_ADDRESS } from 'shared/utils/common/constants';

const axiosInstance = axios.create({
  baseURL: BACKEND_ADDRESS,
});

axiosInstance.interceptors.request.use((config) => {
  const token = window.localStorage.getItem('jwt');

  if (token) {
    // eslint-disable-next-line no-param-reassign
    config.headers.auth = token;
  }
  return config;
});

axiosInstance.interceptors.response.use((response) => response, (error) => {
  if (error?.response?.status === 440) {
    localStorage.removeItem('user');
    localStorage.removeItem('jwt');
    window.location.assign(window.location);
  }

  return Promise.reject(error);
});

export default axiosInstance;
