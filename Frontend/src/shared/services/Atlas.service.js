import axiosInstance from "./axiosInstance"

const MODEL = "atlases"

const AtlasService = {
    getAtlases: async () => {
        const { data } = await axiosInstance.get(`/${MODEL}`)
        return data
    },
}

export default AtlasService