import axiosInstance from './axiosInstance';

const MODEL = 'users';

const UserService = {
  getUsers: async (params) => {
    const { data } = await axiosInstance.get(`/${MODEL}`, { params });
    return data;
  },
};

export default UserService;
