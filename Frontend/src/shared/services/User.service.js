import axiosInstance from './axiosInstance';

const MODEL = 'users';

const UserService = {
  getUsers: async (params) => {
    const preparedParams = { ...params };
    if (!preparedParams.sortBy || preparedParams.sortBy === 'name') {
      preparedParams.sortBy = 'lastName';
    }
    const { data } = await axiosInstance.get(`/${MODEL}`, { params: preparedParams });
    return data;
  },

  getOwnTeams: async () => {
    const { data } = await axiosInstance.get(`/${MODEL}/ownteams`);
    return data;
  },

  getUser: async (id) => {
    const { data } = await axiosInstance.get(`${MODEL}`);
  },
};

export default UserService;
