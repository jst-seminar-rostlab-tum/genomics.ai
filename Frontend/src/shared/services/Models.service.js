import axiosInstance from "./axiosInstance"

const MODEL = "models"

const ModelsService = {
    getModels: async () => {
        const { data } = await axiosInstance.get(`/${MODEL}`)
        return data
    },
}

export default ModelsService