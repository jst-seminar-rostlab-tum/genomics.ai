import React from 'react'
import { Switch, Route, Redirect } from 'react-router-dom'

import Breadcrumb from "components/Breadcrumb"

import { Box } from '@mui/material'
import { TabGroup } from 'components/Tab'
import Search from 'components/Search'
import { useLocation } from 'react-router-dom/cjs/react-router-dom.min'
import Filter from 'components/Filter'

const tmpObj = [
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

  const location = useLocation()

  return (
    <Box>
      {/* Navbar */}
      <Box />

      {/* Breadcrumbs */}
      <Breadcrumb path={`${location.pathname}`} fontSize={1} />

      {/* Search */}
      <Search filterComponent={<Filter references={["test", "test"]} categories={["category1", "category2"]} />} />

      <Box>
        <TabGroup tabsInfo={tmpObj} />
        {/* /explore/atlases */}
        <Switch>
          <Route path="/explore/atlases" render={() => <div>atlases</div>} />
          <Route path="/explore/models" render={() => <div>models</div>} />
          <Route path="/explore/datasets" render={() => <div>datasets</div>} />
          <Redirect to="/explore/atlases" />
        </Switch>
      </Box>
    </Box>
  )
}

export default Explore