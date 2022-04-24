import React from 'react'
import { Switch, Route, Redirect } from 'react-router-dom'

import Breadcrumb from "components/Breadcrumb"

import { Box } from '@mui/material'
import { TabGroup } from 'components/Tab'

const tmpObj=[
  {
    label: "ATLASES",
    path: "/explore/atlases"
  },
  {
    label: "MODELS",
    path: "/explore/models"
  },
  {
    label: "DATASETS",
    path: "/explore/datasets"
  }
]

const Explore = () => {
  return (
    <Box>
        {/* Navbar */}
        <Box />

        {/* Breadcrumbs */}
        <Breadcrumb path="/explore/ajksdflkjasdjkf" fontSize={1} />

        {/* Search */}
        <Box />

        <Box>
            <TabGroup tabsInfo={tmpObj} />
            {/* /explore/atlases */}
            <Switch> 
                <Route path="/explore/atlases" render={()=><div>atlases</div>} />
                <Route path="/explore/models" render={()=><div>models</div>} />
                <Route path="/explore/datasets" render={()=><div>datasets</div>} />
                <Redirect to="/explore/atlases" />
            </Switch>
        </Box>
    </Box>
  )
}

export default Explore