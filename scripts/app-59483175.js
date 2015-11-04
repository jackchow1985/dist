"use strict";angular.module("inspinia",["ngAnimate","ngCookies","ngTouch","ngSanitize","ngResource","ui.router","ui.bootstrap","highcharts-ng","ui.select","angular-venn","angularMoment"]).config(["$stateProvider","$urlRouterProvider",function(t,e){t.state("index",{"abstract":!0,url:"/index",templateUrl:"components/common/content.html"}).state("index.main",{url:"/main",templateUrl:"app/main/main.html",data:{pageTitle:"Example view"}}).state("index.minor",{url:"/minor",templateUrl:"app/minor/minor.html",data:{pageTitle:"Example view"}}),e.otherwise("/index/main")}]),angular.module("inspinia").controller("TabsDemoCtrl",["$scope","$rootScope","$window",function(t){t.tabs=[{title:"Global",content:"Dynamic content 1"},{title:"China",content:"Dynamic content 2",disabled:!0}],t.click=function(){setTimeout(function(){window.dispatchEvent(new Event("resize"))},100)}}]);var app=angular.module("inspinia");app.controller("MapCtrl",["$scope","$location","$http",function(t,e,i){t.data=function(){i.get("getDataFromFile?file=Global_Country_TU").success(function(e){e=e.data.Global_Country_TU;for(var i={},a=0;a<e.length;a++)e[a].Country&&(i[e[a].Country]=parseInt(e[a]["Total Users"]));t.datamap=i})}}]),app.directive("map",function(){return{restrict:"EAC",link:function(t,e){var i=null;t.$watch("datamap",function(){console.log(i),t.datamap&&($(e).width("auto"),$(e).height(400),i=$(e).vectorMap({onRegionClick:function(t,e){console.info(e)},map:"world_mill_en",backgroundColor:"transparent",regionStyle:{initial:{fill:"#e4e4e4","fill-opacity":1,stroke:"none","stroke-width":0,"stroke-opacity":0}},series:{regions:[{values:t.datamap,scale:["#22d6b1","#1ab394"],normalizeFunction:"polynomial"}]},onRegionTipShow:function(e,i,a){i.html(i.html()+": "+(t.datamap[a]?t.datamap[a]:""))}}))})}}}),angular.module("inspinia").controller("VeenCtrl",["$scope","$http","$rootScope",function(t,e,i){t.loadData=function(a){e.get("getDataFromFile?file="+a).success(function(e){e=e.data;var a=[{sets:["China_FX"],size:e.China_FX_FC_SC_TU[0].China_FX},{sets:["China_FC"],size:e.China_FX_FC_SC_TU[0].China_FC},{sets:["China_SC"],size:e.China_FX_FC_SC_TU[0].China_SC},{sets:["China_FX","China_FC"],size:e.China_FX_FC_SC_TU[0].China_FX_FC},{sets:["China_FX","China_SC"],size:e.China_FX_FC_SC_TU[0].China_FX_SC},{sets:["China_FC","China_SC"],size:e.China_FX_FC_SC_TU[0].China_FC_SC},{sets:["China_FX","China_FC","China_SC"],size:e.China_FX_FC_SC_TU[0].China_FX_FC_SC}];i.overlapUsers=e.China_FX_FC_SC_TU[0].China_FX_FC_SC,t.vennData=a,setTimeout(function(){var t=d3.select("#venn"),e=d3.select("body").append("div").attr("class","venntooltip");t.selectAll("g").on("mouseover",function(i){venn.sortAreas(t,i),e.transition().duration(400).style("opacity",.9),e.text(i.size+" users");var a=d3.select(this).transition("tooltip").duration(400);a.select("path").style("stroke-width",3).style("fill-opacity",1==i.sets.length?.4:.1).style("stroke-opacity",1)}).on("mousemove",function(){e.style("left",d3.event.pageX+"px").style("top",d3.event.pageY-28+"px")}).on("mouseout",function(t){e.transition().duration(400).style("opacity",0);var i=d3.select(this).transition("tooltip").duration(400);i.select("path").style("stroke-width",0).style("fill-opacity",1==t.sets.length?.25:0).style("stroke-opacity",0)})},300)})}}]),angular.module("inspinia").controller("SelectorCtrl",["$scope","$rootScope",function(t,e){t.people=[{name:"Adam",email:"adam@email.com",age:12,country:"United States"},{name:"Amalie",email:"amalie@email.com",age:12,country:"Argentina"},{name:"Estefanía",email:"estefania@email.com",age:21,country:"Argentina"},{name:"Adrian",email:"adrian@email.com",age:21,country:"Ecuador"},{name:"Wladimir",email:"wladimir@email.com",age:30,country:"Ecuador"},{name:"Samantha",email:"samantha@email.com",age:30,country:"United States"},{name:"Nicole",email:"nicole@email.com",age:43,country:"Colombia"},{name:"Natasha",email:"natasha@email.com",age:54,country:"Ecuador"},{name:"Michael",email:"michael@email.com",age:15,country:"Colombia"},{name:"Nicolás",email:"nicolas@email.com",age:43,country:"Colombia"}],t.onSelect=function(t){e.$broadcast("changeSelector",t)}}]),angular.module("inspinia").controller("ChartCtrl",["$scope","$http","$rootScope",function(t,e,i){t.$on("changeSelector",function(t,e){console.info(e)}),t.loadData=function(a){e.get("getDataFromFile?file="+a).success(function(e){e=e.data;var n=a.split("|");t.chartConfig.series=[];for(var l=0;l<n.length;l++){for(var o=[],r=0;r<e[n[l]].length;r++)o.push([Date.parse(e[n[l]][r].date),parseInt(e[n[l]][r][t.valueName])]);t.chartConfig.series.push({type:"area",name:t[n[l]],data:o})}if("Date_TU_Global|Date_TU_China"==a){var s=e.Date_TU_Global,c=e.Date_TU_China;i.totalUsers=parseInt(s[s.length-1][t.valueName])+parseInt(c[c.length-1][t.valueName])}console.info(t.chartConfig)})},t.chartConfig={colors:["#A3E1D4","#ff0066","#2b908f","#90ee7e","#eeaaee","#55BF3B","#DF5353","#7798BF","#aaeeee"],chart:{backgroundColor:null,zoomType:"x"},xAxis:{type:"datetime",dateTimeLabelFormats:{month:"%b. %e",year:"%b"},title:{text:"Date"}},yAxis:{title:{text:"Users"},min:0},tooltip:{headerFormat:"<b>{series.name}</b><br>",pointFormat:"{point.x:%e. %b}: {point.y:.f} users"},legend:{enabled:!0},plotOptions:{area:{color:"#0038ff",marker:{enabled:null},lineWidth:1,states:{hover:{lineWidth:1}},threshold:null}},credits:{enabled:!1},series:t.chartSeries}}]),angular.module("inspinia").controller("ChartComLineCtrl",["$scope","$http",function(t,e){t.loadData=function(i){e.get("getDataFromFile?file="+i).success(function(e){e=e.data;var i=[],a=[],n={};n={Others:0};for(var l={},o=0;o<e.Global_Country_TU.length&&4>o;o++)i.push(e.Global_Country_TU[o]),a.push(e.Global_Country_TU[o].Country),n[e.Global_Country_TU[o].Country]=0;for(var o=0;o<e.Date_Global_Country_TU.length;o++)if(a.indexOf(e.Date_Global_Country_TU[o].country)>=0){"TW"==e.Date_Global_Country_TU[o].country&&console.info(e.Date_Global_Country_TU[o]);var r=n[e.Date_Global_Country_TU[o].country]+parseInt(e.Date_Global_Country_TU[o]["Total Users"]);n[e.Date_Global_Country_TU[o].country]=r;var s=[Date.parse(e.Date_Global_Country_TU[o].date),r];l[e.Date_Global_Country_TU[o].country]&&l[e.Date_Global_Country_TU[o].country].length>0?l[e.Date_Global_Country_TU[o].country].push(s):(l[e.Date_Global_Country_TU[o].country]=[],l[e.Date_Global_Country_TU[o].country].push(s))}else{var r=n.Others+parseInt(e.Date_Global_Country_TU[o]["Total Users"]);n.Others=r;var s=[Date.parse(e.Date_Global_Country_TU[o].date),r];l.Others&&l.Others.length>0?l.Others.push(s):(l.Others=[],l.Others.push(s))}var c=[];a.push("Others");for(var h in a)c.push({name:a[h],data:l[a[h]]});t.chartConfig.series=c})},t.chartConfig={chart:{backgroundColor:null,zoomType:"x"},xAxis:{type:"datetime",dateTimeLabelFormats:{month:"%b. %e",year:"%b"},title:{text:"Date"}},yAxis:{title:{text:"Users"},min:0},tooltip:{headerFormat:"<b>{series.name}</b><br>",pointFormat:"{point.x:%e. %b}: {point.y:.f} users"},legend:{enabled:!0},plotOptions:{area:{color:"#0038ff",marker:{enabled:null},lineWidth:1,states:{hover:{lineWidth:1}},threshold:null}},credits:{enabled:!1},series:t.chartSeries}}]),angular.module("inspinia").controller("ChartAreaCtrl",["$scope","$http","$rootScope",function(t,e,i){t.loadData=function(a){e.get("getDataFromFile?file="+a).success(function(e){e=e.data;var n=a.split("|");t.chartConfig.series=[];for(var l=0;l<n.length;l++){for(var o=[],r=0;r<e[n[l]].length;r++)"Date_AU_Global_FX"==a&&i.lastUpdate==e[n[l]][r].date&&(i.activeGlobalFXUsers=parseInt(e[n[l]][r][t.valueName])),"Date_AU_China_FX"==a&&i.lastUpdate==e[n[l]][r].date&&(i.activeChinaFXUsers=parseInt(e[n[l]][r][t.valueName])),"Date_AU_China_FC"==a&&i.lastUpdate==e[n[l]][r].date&&(i.activeChinaFCUsers=parseInt(e[n[l]][r][t.valueName])),"Date_AU_China_SC"==a&&i.lastUpdate==e[n[l]][r].date&&(i.activeChinaSCUsers=parseInt(e[n[l]][r][t.valueName])),"Date_NU_China"==a&&i.lastUpdate==e[n[l]][r].date&&(i.newChinaUsers=parseInt(e[n[l]][r][t.valueName])),"Date_NU_Global"==a&&i.lastUpdate==e[n[l]][r].date&&(i.newGlobalUsers=parseInt(e[n[l]][r][t.valueName])),o.push([Date.parse(e[n[l]][r].date),parseInt(e[n[l]][r][t.valueName])]);t.chartConfig.series.push({type:"area",name:t[n[l]],data:o})}console.info(t.chartConfig)})},t.chartConfig={chart:{backgroundColor:null,zoomType:"x"},xAxis:{type:"datetime",dateTimeLabelFormats:{month:"%b. %e",year:"%b"},title:{text:"Date"}},yAxis:{title:{text:"Users"},min:0},tooltip:{headerFormat:"<b>{series.name}</b><br>",pointFormat:"{point.x:%e. %b}: {point.y:.f} users"},legend:{enabled:!0},plotOptions:{area:{color:"#0038ff",marker:{enabled:null},lineWidth:1,states:{hover:{lineWidth:1}},threshold:null}},credits:{enabled:!1},series:t.chartSeries}}]),angular.module("inspinia").controller("BoxCtrl",["$scope","$http","$rootScope","$filter",function(t,e,i){i.lastUpdate=moment().subtract(2,"days").format("YYYY-MM-DD")}]),$(document).ready(function(){function t(){var t=$("body > #wrapper").height()-61;$(".sidebard-panel").css("min-height",t+"px");var e=$("nav.navbar-default").height(),i=$("#page-wrapper").height();e>i&&$("#page-wrapper").css("min-height",e+"px"),i>e&&$("#page-wrapper").css("min-height",$(window).height()+"px"),$("body").hasClass("fixed-nav")&&$("#page-wrapper").css("min-height",$(window).height()-60+"px")}$(window).bind("load resize scroll",function(){$("body").hasClass("body-small")||t()}),setTimeout(function(){t()})}),$(function(){$(window).bind("load resize",function(){$(this).width()<769?$("body").addClass("body-small"):$("body").removeClass("body-small")})}),angular.module("inspinia").directive("sideNavigation",["$timeout",function(t){return{restrict:"A",link:function(e,i){e.$watch("authentication.user",function(){t(function(){i.metisMenu()})})}}}]).directive("minimalizaSidebar",["$timeout",function(t){return{restrict:"A",template:'<a class="navbar-minimalize minimalize-styl-2 btn btn-primary " href="" ng-init="minimalize()"><i class="fa fa-bars"></i></a>',controller:["$scope","$element",function(e){e.minimalize=function(){angular.element("body").toggleClass("mini-navbar"),!angular.element("body").hasClass("mini-navbar")||angular.element("body").hasClass("body-small")?(angular.element("#side-menu").hide(),t(function(){angular.element("#side-menu").fadeIn(500)},100)):angular.element("#side-menu").removeAttr("style")}}]}}]),angular.module("inspinia").run(["$templateCache",function(t){t.put("app/main/main.html",'<div class="wrapper wrapper-content animated fadeInRight"><div ng-include="\'components/box/sumbox.html\'"></div><div ng-include="\'components/chart/chart_tu.html\'"></div><div class="row chartPanel"><div ng-include="\'components/map/map.html\'"></div><div ng-include="\'components/chart/chart_global_dist.html\'"></div></div><div class="row chartPanel"><div ng-include="\'components/chart/chart_fx_fc_sc.html\'"></div></div><div class="row chartPanel"><div ng-include="\'components/chart/overlap.html\'"></div><div ng-include="\'components/chart/veen.html\'"></div></div><div ng-controller="TabsDemoCtrl" class="tab-bar"><tabset><tab heading="Global"><div class="row chartPanel"><div ng-include="\'components/chart/chart_gb_nu.html\'"></div></div></tab><tab heading="China" ng-click="click()"><div class="row chartPanel"><div ng-include="\'components/chart/chart_cn_nu.html\'"></div></div></tab></tabset></div><div ng-controller="TabsDemoCtrl" class="tab-bar"><tabset><tab heading="Global"><div class="row chartPanel"><div ng-include="\'components/chart/chart_gb_au.html\'"></div></div></tab><tab heading="China FX" ng-click="click()"><div class="row chartPanel"><div ng-include="\'components/chart/chart_cn_au_fx.html\'"></div></div></tab><tab heading="China FC" ng-click="click()"><div class="row chartPanel"><div ng-include="\'components/chart/chart_cn_au_fc.html\'"></div></div></tab><tab heading="China SC" ng-click="click()"><div class="row chartPanel"><div ng-include="\'components/chart/chart_cn_au_sc.html\'"></div></div></tab></tabset></div></div>'),t.put("app/minor/minor.html",'<div class="wrapper wrapper-content animated fadeInRight"><div class="row"><div class="col-lg-12"><div class="text-center m-t-lg"><h1>Simple example of second view</h1><small>Configure in app.js as index.minor state.</small></div></div></div></div>'),t.put("components/box/sumbox.html",'<div class="row" ng-controller="BoxCtrl"><div class="col-lg-3"><div class="ibox float-e-margins"><div class="ibox-title"><span class="label label-success pull-right">All Apps</span><h5>Total Users</h5></div><div class="ibox-content"><h1 class="no-margins" ng-bind="($root.totalUsers|number)">0</h1><small>Users</small></div></div></div><div class="col-lg-3"><div class="ibox float-e-margins"><div class="ibox-title"><span class="label label-info pull-right">Daily</span><h5>Active User</h5></div><div class="ibox-content"><h1 class="no-margins" ng-bind="(($root.activeChinaFXUsers + $root.activeChinaSCUsers + $root.activeChinaFCUsers + $root.activeGlobalFXUsers ) |number)">0</h1><div class="stat-percent font-bold text-success" ng-bind="(($root.activeChinaFXUsers + $root.activeChinaSCUsers + $root.activeChinaFCUsers + $root.activeGlobalFXUsers ) / $root.totalUsers * 100|number:2) + \'% active\'"><i class="fa fa-bolt"></i></div><small>Users</small></div></div></div><div class="col-lg-3"><div class="ibox float-e-margins"><div class="ibox-title"><span class="label label-primary pull-right">Daily</span><h5>New Users</h5></div><div class="ibox-content"><h1 class="no-margins" ng-bind="($root.overlapUsers |number)">0</h1><div class="stat-percent font-bold text-info" ng-bind="(($root.overlapUsers ) / $root.totalUsers * 100|number:2) + \'% new\'"><i class="fa fa-bolt"></i></div><small>Users</small></div></div></div><div class="col-lg-3"><div class="ibox float-e-margins"><div class="ibox-title"><span class="label label-danger pull-right">China</span><h5>Overlap Users</h5></div><div class="ibox-content"><h1 class="no-margins">80,600</h1><small>Users</small></div></div></div></div>'),t.put("components/common/content.html",'<div id="wrapper"><div ng-include="\'components/common/navigation.html\'"></div><div id="page-wrapper" class="gray-bg {{$state.current.name}}"><div ng-include="\'components/common/topnavbar.html\'"></div><div ui-view=""></div><div ng-include="\'components/common/footer.html\'"></div></div></div>'),t.put("components/common/footer.html",'<div class="footer"><div class="pull-right">10GB of <strong>250GB</strong> Free.</div><div><strong>Copyright</strong> Example Company &copy; 2014-2015</div></div>'),t.put("components/common/ibox_tools.html",'<div class="ibox-tools dropdown" dropdown=""><a ng-click="showhide()"><i class="fa fa-chevron-up"></i></a> <a class="dropdown-toggle" href="" dropdown-toggle=""><i class="fa fa-wrench"></i></a><ul class="dropdown-menu dropdown-user"><li><a href="">Config option 1</a></li><li><a href="">Config option 2</a></li></ul><a ng-click="closebox()"><i class="fa fa-times"></i></a></div>'),t.put("components/common/navigation.html",'<nav class="navbar-default navbar-static-side" role="navigation"><div class="sidebar-collapse"><ul side-navigation="" class="nav metismenu" id="side-menu"><li class="nav-header"><div class="dropdown profile-element" dropdown=""><a class="dropdown-toggle" dropdown-toggle="" href=""><span class="clear"><span class="block m-t-xs"><strong class="font-bold">{{main.userName}}</strong></span> <span class="text-muted text-xs block">Example menu<b class="caret"></b></span></span></a><ul class="dropdown-menu animated fadeInRight m-t-xs"><li><a href="">Logout</a></li></ul></div><div class="logo-element">FDT</div></li><li ui-sref-active="active"><a ui-sref="index.main"><i class="fa fa-users"></i> <span class="nav-label">User</span></a></li><li ui-sref-active="active"><a ui-sref="index.school"><i class="fa fa-graduation-cap"></i> <span class="nav-label">School</span></a></li><li ui-sref-active="active"><a ui-sref="index.region"><i class="fa fa-map-marker"></i> <span class="nav-label">Region</span></a></li></ul></div></nav>'),t.put("components/common/topnavbar.html",'<div class="row border-bottom"><nav class="navbar navbar-static-top white-bg" role="navigation" style="margin-bottom: 0"><div class="navbar-header"><span minimaliza-sidebar=""></span><form role="search" class="navbar-form-custom" method="post" action=""><div class="form-group"><input type="text" placeholder="Search for something..." class="form-control" name="top-search" id="top-search"></div></form></div><ul class="nav navbar-top-links navbar-right"><li><span class="m-r-sm text-muted welcome-message" ng-bind="\'Last Updated Time : \' + $root.lastUpdate"></span></li><li><a href=""><i class="fa fa-sign-out"></i> Log out</a></li></ul></nav></div>'),t.put("components/map/map.html",'<div class="col-md-6 chartPanel" ng-controller="MapCtrl"><h3>Global Distribution on Map</h3><div id="map" map="" ng-init="data()"></div></div>'),t.put("components/chart/chart_cn_au_fc.html",'<div class="col-md-4" ng-controller="ChartAreaCtrl"><highchart ng-init="loadData(\'Date_AU_China_FC\')" config="chartConfig"></highchart><input type="hidden" ng-init="chartConfig.title.text=\'Daily Active Users\'"> <input type="hidden" ng-init="chartConfig.subtitle.text=\'China\'"> <input type="hidden" ng-init="valueName=\'Active Users\'"> <input type="hidden" ng-init="Date_AU_China_FC=\'China Daily Active Users\'"></div><div class="col-md-4" ng-controller="ChartAreaCtrl"><highchart ng-init="loadData(\'Week_AU_China_FC\')" config="chartConfig"></highchart><input type="hidden" ng-init="chartConfig.title.text=\'Weekly Active Users\'"> <input type="hidden" ng-init="chartConfig.subtitle.text=\'China\'"> <input type="hidden" ng-init="valueName=\'Active Users\'"> <input type="hidden" ng-init="Week_AU_China_FC=\'China Weekly Active Users\'"></div><div class="col-md-4" ng-controller="ChartAreaCtrl"><highchart ng-init="loadData(\'Month_AU_China_FC\')" config="chartConfig"></highchart><input type="hidden" ng-init="chartConfig.title.text=\'Monthly Active Users\'"> <input type="hidden" ng-init="chartConfig.subtitle.text=\'China\'"> <input type="hidden" ng-init="valueName=\'Active Users\'"> <input type="hidden" ng-init="Month_AU_China_FC=\'China Monthly Active Users\'"></div>'),t.put("components/chart/chart_cn_au_fx.html",'<div class="col-md-4" ng-controller="ChartAreaCtrl"><highchart ng-init="loadData(\'Date_AU_China_FX\')" config="chartConfig"></highchart><input type="hidden" ng-init="chartConfig.title.text=\'Daily Active Users\'"> <input type="hidden" ng-init="chartConfig.subtitle.text=\'China\'"> <input type="hidden" ng-init="valueName=\'Active Users\'"> <input type="hidden" ng-init="Date_AU_China_FX=\'China Daily Active Users\'"></div><div class="col-md-4" ng-controller="ChartAreaCtrl"><highchart ng-init="loadData(\'Week_AU_China_FX\')" config="chartConfig"></highchart><input type="hidden" ng-init="chartConfig.title.text=\'Weekly Active Users\'"> <input type="hidden" ng-init="chartConfig.subtitle.text=\'China\'"> <input type="hidden" ng-init="valueName=\'Active Users\'"> <input type="hidden" ng-init="Week_AU_China_FX=\'China Weekly Active Users\'"></div><div class="col-md-4" ng-controller="ChartAreaCtrl"><highchart ng-init="loadData(\'Month_AU_China_FX\')" config="chartConfig"></highchart><input type="hidden" ng-init="chartConfig.title.text=\'Monthly Active Users\'"> <input type="hidden" ng-init="chartConfig.subtitle.text=\'China\'"> <input type="hidden" ng-init="valueName=\'Active Users\'"> <input type="hidden" ng-init="Month_AU_China_FX=\'China Monthly Active Users\'"></div>'),t.put("components/chart/chart_cn_au_sc.html",'<div class="col-md-4" ng-controller="ChartAreaCtrl"><highchart ng-init="loadData(\'Date_AU_China_SC\')" config="chartConfig"></highchart><input type="hidden" ng-init="chartConfig.title.text=\'Daily Active Users\'"> <input type="hidden" ng-init="chartConfig.subtitle.text=\'China\'"> <input type="hidden" ng-init="valueName=\'Active Users\'"> <input type="hidden" ng-init="Date_AU_China_SC=\'China Daily Active Users\'"></div><div class="col-md-4" ng-controller="ChartAreaCtrl"><highchart ng-init="loadData(\'Week_AU_China_SC\')" config="chartConfig"></highchart><input type="hidden" ng-init="chartConfig.title.text=\'Weekly Active Users\'"> <input type="hidden" ng-init="chartConfig.subtitle.text=\'China\'"> <input type="hidden" ng-init="valueName=\'Active Users\'"> <input type="hidden" ng-init="Week_AU_China_SC=\'China Weekly Active Users\'"></div><div class="col-md-4" ng-controller="ChartAreaCtrl"><highchart ng-init="loadData(\'Month_AU_China_SC\')" config="chartConfig"></highchart><input type="hidden" ng-init="chartConfig.title.text=\'Monthly Active Users\'"> <input type="hidden" ng-init="chartConfig.subtitle.text=\'China\'"> <input type="hidden" ng-init="valueName=\'Active Users\'"> <input type="hidden" ng-init="Month_AU_China_SC=\'China Monthly Active Users\'"></div>'),t.put("components/chart/chart_cn_nu.html",'<div class="col-md-4" ng-controller="ChartAreaCtrl"><highchart ng-init="loadData(\'Date_NU_China\')" config="chartConfig"></highchart><input type="hidden" ng-init="chartConfig.title.text=\'Daily New Users\'"> <input type="hidden" ng-init="valueName=\'New Users\'"> <input type="hidden" ng-init="Date_NU_China=\'China Daily New Users\'"></div><div class="col-md-4" ng-controller="ChartAreaCtrl"><highchart ng-init="loadData(\'Week_NU_China\')" config="chartConfig"></highchart><input type="hidden" ng-init="chartConfig.title.text=\'Weekly New Users\'"> <input type="hidden" ng-init="valueName=\'New Users\'"> <input type="hidden" ng-init="Week_NU_China=\'China Weekly New Users\'"></div><div class="col-md-4" ng-controller="ChartAreaCtrl"><highchart ng-init="loadData(\'Month_NU_China\')" config="chartConfig"></highchart><input type="hidden" ng-init="chartConfig.title.text=\'Monthly New Users\'"> <input type="hidden" ng-init="valueName=\'New Users\'"> <input type="hidden" ng-init="Month_NU_China=\'China Monthly New Users\'"></div>'),t.put("components/chart/chart_fx_fc_sc.html",'<div class="col-md-12" ng-controller="ChartCtrl"><highchart ng-init="loadData(\'Date_TU_China_FX|Date_TU_China_FC|Date_TU_China_SC\')" config="chartConfig"></highchart><input type="hidden" ng-init="chartConfig.title.text=\'China FX FC SC\'"> <input type="hidden" ng-init="valueName=\'Total Users\'"> <input type="hidden" ng-init="Date_TU_China_FX=\'China FX\'"> <input type="hidden" ng-init="Date_TU_China_FC=\'China FC\'"> <input type="hidden" ng-init="Date_TU_China_SC=\'China SC\'"></div>'),t.put("components/chart/chart_gb_au.html",'<div class="col-md-4" ng-controller="ChartAreaCtrl"><highchart ng-init="loadData(\'Date_AU_Global_FX\')" config="chartConfig"></highchart><input type="hidden" ng-init="chartConfig.title.text=\'Daily Active Users\'"> <input type="hidden" ng-init="chartConfig.subtitle.text=\'Global\'"> <input type="hidden" ng-init="valueName=\'Active Users\'"> <input type="hidden" ng-init="Date_AU_Global_FX=\'Global Daily Active Users\'"></div><div class="col-md-4" ng-controller="ChartAreaCtrl"><highchart ng-init="loadData(\'Week_AU_Global_FX\')" config="chartConfig"></highchart><input type="hidden" ng-init="chartConfig.title.text=\'Weekly Active Users\'"> <input type="hidden" ng-init="chartConfig.subtitle.text=\'Global\'"> <input type="hidden" ng-init="valueName=\'Active Users\'"> <input type="hidden" ng-init="Week_AU_Global_FX=\'Global Weekly Active Users\'"></div><div class="col-md-4" ng-controller="ChartAreaCtrl"><highchart ng-init="loadData(\'Month_AU_Global_FX\')" config="chartConfig"></highchart><input type="hidden" ng-init="chartConfig.title.text=\'Monthly Active Users\'"> <input type="hidden" ng-init="chartConfig.subtitle.text=\'Global\'"> <input type="hidden" ng-init="valueName=\'Active Users\'"> <input type="hidden" ng-init="Month_AU_Global_FX=\'Global Monthly Active Users\'"></div>'),t.put("components/chart/chart_gb_nu.html",'<div class="col-md-4" ng-controller="ChartAreaCtrl"><highchart ng-init="loadData(\'Date_NU_Global\')" config="chartConfig"></highchart><input type="hidden" ng-init="chartConfig.title.text=\'Daily New Users\'"> <input type="hidden" ng-init="chartConfig.subtitle.text=\'Global\'"> <input type="hidden" ng-init="valueName=\'New Users\'"> <input type="hidden" ng-init="Date_NU_Global=\'Global Daily New Users\'"></div><div class="col-md-4" ng-controller="ChartAreaCtrl"><highchart ng-init="loadData(\'Week_NU_Global\')" config="chartConfig"></highchart><input type="hidden" ng-init="chartConfig.title.text=\'Weekly New Users\'"> <input type="hidden" ng-init="chartConfig.subtitle.text=\'Global\'"> <input type="hidden" ng-init="valueName=\'New Users\'"> <input type="hidden" ng-init="Week_NU_Global=\'Global Weekly New Users\'"></div><div class="col-md-4" ng-controller="ChartAreaCtrl"><highchart ng-init="loadData(\'Month_NU_Global\')" config="chartConfig"></highchart><input type="hidden" ng-init="chartConfig.title.text=\'Monthly New Users\'"> <input type="hidden" ng-init="chartConfig.subtitle.text=\'Global\'"> <input type="hidden" ng-init="valueName=\'New Users\'"> <input type="hidden" ng-init="Month_NU_Global=\'Global Monthly New Users\'"></div>'),t.put("components/chart/chart_global_dist.html",'<div ng-controller="ChartComLineCtrl"><div class="col-md-6"><highchart ng-init="loadData(\'Date_Global_Country_TU|Global_Country_TU\')" config="chartConfig"></highchart><input type="hidden" ng-init="chartConfig.title.text=\'Global Distribution\'"></div></div>'),t.put("components/chart/chart_tu.html",'<div class="row" ng-controller="SelectorCtrl"><div class="col-md-3 selector pull-right"><ui-select ng-model="person.selected" on-select="onSelect($item, $model)" theme="select2" class="form-control" title="Choose a person"><ui-select-match placeholder="Select or search a person in the list...">{{$select.selected.name}}</ui-select-match><ui-select-choices repeat="item in people | filter: $select.search"><div ng-bind-html="item.name | highlight: $select.search"></div><small ng-bind-html="item.email | highlight: $select.search"></small></ui-select-choices></ui-select></div></div><div class="row" ng-controller="ChartCtrl"><div class="col-lg-12"><highchart ng-init="loadData(\'Date_TU_Global|Date_TU_China\')" config="chartConfig" class="span9"></highchart><input type="hidden" ng-init="chartConfig.title.text=\'Total User Growth\'"> <input type="hidden" ng-init="valueName=\'Total Users\'"> <input type="hidden" ng-init="Date_TU_Global=\'Global\'"> <input type="hidden" ng-init="Date_TU_China=\'China\'"> <input type="hidden" ng-init="chartType=\'area\'"> <input type="hidden" ng-init="chartConfig.colors=[\'#0000ff\',\'ff0066\',\'2b908f\']"></div></div>'),t.put("components/chart/overlap.html",'<div class="col-md-6" ng-controller="ChartCtrl"><highchart ng-init="loadData(\'FX_FC_Overlap|FX_SC_Overlap|FC_SC_Overlap|FX_FC_SC_Overlap\')" config="chartConfig"></highchart><input type="hidden" ng-init="chartConfig.title.text=\'FX FC SC Users\'"> <input type="hidden" ng-init="chartConfig.subtitle.text=\'Overlaping Users\'"> <input type="hidden" ng-init="valueName=\'Common Users\'"> <input type="hidden" ng-init="FX_FC_Overlap=\'FX&FC\'"> <input type="hidden" ng-init="FX_SC_Overlap=\'FX&SC\'"> <input type="hidden" ng-init="FC_SC_Overlap=\'FC&SC\'"> <input type="hidden" ng-init="FX_FC_SC_Overlap=\'FX&FC&SC\'"></div>'),t.put("components/chart/veen.html",'<div class="col-md-6 venn-panel" ng-controller="VeenCtrl"><div id="venn" venn="vennData" ng-init="loadData(\'China_FX_FC_SC_TU\')"></div></div>')}]);