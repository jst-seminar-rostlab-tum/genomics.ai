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

export default axiosInstance;
