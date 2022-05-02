import axiosInstance from "./axiosInstance"

const MODEL = "atlases"

const AtlasService = {
  getAtlases: async () => {
    const { data } = await axiosInstance.get(`/${MODEL}`)
    return data
  },
  getAtlasById: async (id) => {
    const { data } = await axiosInstance.get(`/${MODEL}/${id}`)
    return data
  }
}

export default AtlasService