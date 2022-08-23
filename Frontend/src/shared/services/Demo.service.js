import axiosInstance from './axiosInstance';

const DEMOS = 'demos';

const DemoService = {
  getDemos: async () => {
    const { data } = await axiosInstance.get(`/${DEMOS}`);
    return data;
  },
};

export default DemoService;
