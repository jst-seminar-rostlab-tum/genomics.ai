import React, { useEffect } from 'react'
import AtlasService from 'shared/services/Atlas.service'
import ModelsService from 'shared/services/Models.service'

const MapperContext = React.createContext()
MapperContext.displayName = 'MapperContext'

const MapperProvider = (props) => {
  const [atlasId, setAtlasId] = useState(localStorage.getItem("atlasId") ? JSON.parse(localStorage.getItem("atlasId")) : null) 
  const [modelId, setModelId] = useState(localStorage.getItem("modelId") ? JSON.parse(localStorage.getItem("modelId")) : null) 
  const [selectedAtlas, setSelectedAtlas] = useState(null)
  const [selectedModel, setSelectedModel] = useState(null)

  useEffect(() => {
    if(atlasId) {
      AtlasService.getAtlasById(atlasId)
        .then((data) => setSelectedAtlas(data))
        .catch((err) => console.log(err))
    }
    if(modelId) {
      ModelsService.getModelById(modelId)
        .then((data) => setSelectedModel(data))
        .catch((err) => console.log(err))
    }
  }, [atlasId, modelId])

  const value = [{setAtlasId, setModelId}, {selectedAtlas, selectedModel}] 

  return <MapperContext.Provider value={value} {...props}/>
}
